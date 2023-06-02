import numpy as np
import imageProcessing as imProc
from matplotlib.widgets import Slider, Button, RadioButtons

class interactivePlot(object):   
    def __init__(self,fig, ax, X, fil, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
#        set various event handlers
        fig.canvas.mpl_connect('scroll_event', self.onScroll)
        fig.canvas.mpl_connect('button_press_event', self.onClick)
        fig.canvas.mpl_connect('button_release_event', self.onRelease)
        fig.canvas.mpl_connect('motion_notify_event', self.onMotion)
        fig.canvas.mpl_connect('key_press_event', self.keyPress)
        
        self.fig = fig
        self.plotAxis = plotAxis
        self.ax = ax
        self.ax.set_adjustable('box')
        self.mouseClicked = None
        self.cmapRange = np.nanmax(X) - np.nanmin(X)
        self.cmapCenter = np.nanmin(X) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.X = X
        self.originallog = np.log(self.X)
        self.axisLabels = axisLabels
        self.colorMap = cmap
        self.filter_val = 0
        self.original = X
        self.fil = fil
        
        if fov is None:
            self.fov = (1,1,1)
            self.resolution = (1,1,1)
        else:
            self.fov = fov
            self.resolution = np.divide(fov, np.shape(X))

        self.slices = np.shape(X)[plotAxis]
        self.ind = self.slices//2
        
        ax.set_title('Slice %s/%s' % (self.ind,self.slices))
        
        if self.plotAxis == 0:
            imageData = self.X[self.ind,:,:]
        elif self.plotAxis == 1:
            imageData = self.X[:,self.ind,:]
        elif self.plotAxis == 2:
            imageData = self.X[:,:,self.ind]
        else:
            print("invalid axis")
            return -1
        self.im = ax.imshow(imageData, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        
        #create slider
        self.axSlider = fig.add_axes([0.15, 0.03, 0.6, 0.03])
        self.slider = Slider(ax=self.axSlider, label='', valmin=0, valmax=1, valinit=0)
        
        def update(sliderVal):
            if self.log == 'raw data':
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(self.original,sliderVal)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(self.original,sliderVal)
                self.X = kSpace
                self.updateSlice()
            elif self.log == 'log(data)':
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(self.originallog,sliderVal)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(self.originallog,sliderVal)
                self.X = kSpace
                self.updateSlice()

        self.slider.on_changed(update)
        self.updateSlice()
        
        #create radiobutton
        self.axRadio = fig.add_axes([0.85, 0.8, 0.15, 0.1])
        self.radio = RadioButtons(self.axRadio, ('raw data', 'log(data)'))
        self.log = 'raw data'
    
        def data_transform(label):
            self.log = label
            if label == 'log(data)':
                self.X = np.log(self.X)
            elif label == 'raw data':
                self.X = np.exp(self.X)
            self.updateSlice()

        self.radio.on_clicked(data_transform)

    def keyPress(self, event):
        if event.key == " ": #change orientation on space bar press
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
        elif event.key == "up": #change orientation on space bar press
            self.ind = (self.ind + 1) % self.slices
        elif event.key == "down": #change orientation on space bar press
            self.ind = (self.ind - 1) % self.slices
        elif event.key == "right":
            if self.filter_val < 1:
                self.filter_val += 1/100
                kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.original))))
                
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(kSpace,self.filter_val)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(kSpace,self.filter_val)
                    
                self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))
        elif event.key == "left":
            if self.filter_val > 0.01:
                self.filter_val -= 1/100
                kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.original))))
                
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(kSpace,self.filter_val)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(kSpace,self.filter_val)
                
                self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))    
        self.updateSlice()


    def onScroll(self, event):
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.updateSlice()
    
    def onClick(self, event):
#        print('Button pressed: %s at X: %s, Y: %s'%(event.button, event.xdata, event.ydata))
        if event.button == 3:   #Reset image on right click
            self.cmapRange = np.nanmax(self.X) - np.nanmin(self.X)
            self.cmapCenter = np.nanmin(self.X) + self.cmapRange/2
            self.tempCmapRange = self.cmapRange
            self.tempCmapCenter = self.cmapCenter
            self.im.set_clim(self.cmapCenter-self.cmapRange/2, self.cmapCenter+self.cmapRange/2)
            self.im.axes.figure.canvas.draw()
        elif event.button == 2: #change orientation on scroll wheel click
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.updateSlice()
        else:
            self.mouseClicked = event.xdata, event.ydata
        
    def onRelease(self, event):
        self.mouseClicked = None
        self.cmapRange = self.tempCmapRange
        self.cmapCenter = self.tempCmapCenter
        
    def onMotion(self, event):
        if self.mouseClicked == None: return        #if mouse isn't clicked ignore movement

        dx = event.xdata - self.mouseClicked[0]
        dy = event.ydata - self.mouseClicked[1]
        
        normDx = dx/self.mouseClicked[0]
        normDy = dy/self.mouseClicked[1]
        
        self.tempCmapRange = self.cmapRange*(1+normDy)
        self.tempCmapCenter = self.cmapCenter*(1+normDx)
        
        self.im.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
        self.im.axes.figure.canvas.draw()
        
    def updateSlice(self):
        if self.plotAxis == 0:
            if self.ind >= np.size(self.X,0):
                self.ind = np.size(self.X, 0)-1
            imageData = self.X[self.ind,:,:]
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[1])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.ax.set_aspect(self.resolution[1]/self.resolution[2])
            self.slices = np.size(self.X,0)
        elif self.plotAxis == 1:
            if self.ind >= np.size(self.X,1):
                self.ind = np.size(self.X, 1)-1
            imageData = self.X[:,self.ind,:]
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[0])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.ax.set_aspect(self.resolution[0]/self.resolution[2])
            self.slices = np.size(self.X,1)
        else:
            if self.ind >= np.size(self.X,2):
                self.ind = np.size(self.X,2)-1
            imageData = self.X[:,:,self.ind]
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[1])
                self.ax.set_ylabel(self.axisLabels[0])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.ax.set_aspect(self.resolution[0]/self.resolution[1])
            self.slices = np.size(self.X,2)
        # print("Plot axis: %.0f, Array dims: %.0f , %.0f" %(self.plotAxis, np.size(imageData,0), np.size(imageData,1)))
        #self.im.set_data(imageData)
        self.cmapRange = np.nanmax(self.X) - np.nanmin(self.X)
        self.cmapCenter = np.nanmin(self.X) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.im = self.ax.imshow(imageData, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax.set_title('Slice %s/%s' % (self.ind,self.slices))
        self.im.axes.figure.canvas.draw()
        
