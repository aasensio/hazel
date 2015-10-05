"""
This demo demonstrates how to embed a matplotlib (mpl) plot 
into a PyQt4 GUI application, including:

* Using the navigation toolbar
* Adding data to the plot
* Dynamically modifying the plot's properties
* Processing mpl events
* Saving the plot to a file from a menu

The main goal is to serve as a basis for developing rich PyQt GUI
applications featuring mpl plots (using the mpl OO API).

Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 19.01.2009
"""
import sys, os, random
import numpy as np
import i0Allen
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

class FloatSlider(QSlider):
    def __init__(self, parent, decimals=3, *args, **kargs):
        super(FloatSlider, self).__init__(parent, *args, **kargs)
        self._multi = 10 ** decimals
        self.setMinimum(self.minimum())
        self.setMaximum(self.maximum())

    def value(self):
        return float(super(FloatSlider, self).value()) / self._multi

    def setMinimum(self, value):
        return super(FloatSlider, self).setMinimum(value * self._multi)

    def setMaximum(self, value):
        return super(FloatSlider, self).setMaximum(value * self._multi)

    def setValue(self, value):
        super(FloatSlider, self).setValue(int(value * self._multi))

class ExtendedQLabel(QLabel):
 
    def __init(self, parent):
        QLabel.__init__(self, parent)
 
    def mouseReleaseEvent(self, ev):
        self.emit(SIGNAL('clicked()'))


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('hazel')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        
        # self.on_draw()
    
#######################################################################
# EVENTS
#######################################################################
    def on_about(self):
        msg = """ Hazel:
        
         * A. Asensio Ramos
         * Instituto de Astrofisica de Canarias, Spain
        """
        QMessageBox.about(self, "About hazel", msg.strip())
    
    # def on_pick(self, event):
    #     # The event received here is of the type
    #     # matplotlib.backend_bases.PickEvent
    #     #
    #     # It carries lots of information, of which we're using
    #     # only a small amount here.
    #     # 
    #     box_points = event.artist.get_bbox().get_points()
    #     msg = "You've clicked on a bar with coords:\n %s" % box_points
        
    #     QMessageBox.information(self, "Click!", msg)
    
    # def on_draw(self):
    #     """ Redraws the figure
    #     """
    #     self.data = np.arange(10)
        
    #     x = np.arange(10)

    #     # clear the axes and redraw the plot anew
    #     #
    #     self.axes.clear()        
    #     # self.axes.grid(self.grid_cb.isChecked())
        
    #     self.axes.bar(
    #         left=x, 
    #         height=self.data, 
    #         width=self.slider.value() / 100.0, 
    #         align='center', 
    #         alpha=0.44,
    #         picker=5)
        
    #     self.canvas.draw()

    def onSliderB1(self):
        self.sliderValueB1.setText('{0:6.1f}'.format(self.sliderB1.value()))

    def onSliderthB1(self):
        self.sliderValuethB1.setText(str(self.sliderthB1.value()))

    def onSliderLeftB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter lower limit:')
        if (ok):
            self.sliderB1.setMinimum(float(text))
            self.sliderLeftB1.setText(text)

    def onSliderRightB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter upper limit:')
        if (ok):
            self.sliderB1.setMaximum(float(text))
            self.sliderRightB1.setText(text)

    def onSliderValueB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value:')
        if (ok):
            self.sliderB1.setValue(float(text))
            self.sliderValueB1.setText(text)

    def onSliderB2(self):
        self.sliderValueB2.setText(str(self.sliderB2.value()))

    def onSliderthB2(self):
        self.sliderValuethB2.setText(str(self.sliderthB2.value()))

    def onSliderLeftB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter lower limit:')
        if (ok):
            self.sliderB2.setMinimum(float(text))
            self.sliderLeftB2.setText(text)

    def onSliderRightB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter upper limit:')
        if (ok):
            self.sliderB2.setMaximum(float(text))
            self.sliderRightB2.setText(text)

    def onSliderTheta(self):
        self.sliderValuertheta.setText(str(self.slidertheta.value()))
        if (self.checkAllen.isChecked()):
            theta = float(self.sliderValuertheta.text())
            mu = np.cos(theta * np.pi / 180.0)
            i0 = i0Allen.i0Allen(10830.0, mu)
            self.I0.setText(str(i0))

    def onSliderValueTheta(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter theta angle:')
        if (ok):
            self.slidertheta.setValue(float(text))
            self.sliderValuertheta.setText(str(text))

    def onSliderPhi(self):
        self.sliderValuerphi.setText(str(self.sliderphi.value()))

    def onSliderValuePhi(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter phi angle:')
        if (ok):
            self.sliderphi.setValue(float(text))
            self.sliderValuerphi.setText(str(text))

    def onSliderGamma(self):
        self.sliderValuergamma.setText(str(self.slidergamma.value()))

    def onSliderValueGamma(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter gamma angle:')
        if (ok):
            self.slidergamma.setValue(float(text))
            self.sliderValuergamma.setText(str(text))

    def onSliderHeight(self):
        self.sliderValuerheight.setText(str(self.sliderheight.value()))

    def onSliderValueHeight(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter height:')
        if (ok):
            self.sliderheight.setValue(float(text))
            self.sliderValuerheight.setText(str(text))

    def onSliderDamp(self):
        self.sliderValuerDamp.setText(str(self.sliderDamp.value()))

    def onSliderValueDamp(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter damping:')
        if (ok):
            self.sliderDamp.setValue(float(text))
            self.sliderValuerDamp.setText(str(text))

    def onSliderff(self):
        self.sliderValuerff.setText(str(self.sliderff.value()))

    def onSliderValueff(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter ff:')
        if (ok):
            self.sliderff.setValue(float(text))
            self.sliderValuerff.setText(str(text))

    def onSliderS2S1(self):
        self.sliderValuerS2S1.setText(str(self.sliderS2S1.value()))

    def onSliderValueS2S1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter S2/S1:')
        if (ok):
            self.sliderS2S1.setValue(float(text))
            self.sliderValuerS2S1.setText(str(text))

    def onCheckAllen(self):
        self.onSliderTheta()


    def onChangeSlabs(self):
        pass
    
#######################################################################
# INITIALIZATION
#######################################################################
    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
                
        self.axes = self.fig.add_subplot(2,2,1)
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        # self.canvas.mpl_connect('pick_event', self.on_pick)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # First component
        # 
        boxB1 = QGridLayout()        

        sliderLabel = QLabel('B [G]')
        self.sliderB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftB1 = ExtendedQLabel('0')
        self.sliderRightB1 = ExtendedQLabel('1000')
        self.sliderB1.setMinimum(0.0)
        self.sliderB1.setMaximum(1000.0)
        self.sliderB1.setValue(20)
        self.sliderB1.setTracking(True)
        self.sliderB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderValueB1 = ExtendedQLabel('{0:6.1f}'.format(0))        
        for l, w in enumerate([sliderLabel, self.sliderLeftB1, self.sliderB1, self.sliderRightB1, self.sliderValueB1]):
            boxB1.addWidget(w, 0, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderB1, SIGNAL('sliderReleased()'), self.onSliderB1)
        self.connect(self.sliderLeftB1, SIGNAL('clicked()'), self.onSliderLeftB1)
        self.connect(self.sliderRightB1, SIGNAL('clicked()'), self.onSliderRightB1)
        self.connect(self.sliderValueB1, SIGNAL('clicked()'), self.onSliderValueB1)
                

        sliderLabel = QLabel('thB [deg]')
        self.sliderthB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftthB1 = ExtendedQLabel('0')
        self.sliderRightthB1 = ExtendedQLabel('180')
        self.sliderthB1.setMinimum(0.0)
        self.sliderthB1.setMaximum(180.0)
        self.sliderthB1.setValue(20)
        self.sliderthB1.setTracking(True)
        self.sliderthB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuethB1 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderLeftthB1, self.sliderthB1, self.sliderRightthB1, self.sliderValuethB1]):
            boxB1.addWidget(w, 1, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderthB1, SIGNAL('sliderReleased()'), self.onSliderthB1)

        sliderLabel = QLabel('phiB [deg]')
        self.sliderphiB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftphiB1 = ExtendedQLabel('-180')
        self.sliderRightphiB1 = ExtendedQLabel('180')
        self.sliderphiB1.setMinimum(-180)
        self.sliderphiB1.setMaximum(180)
        self.sliderphiB1.setValue(20)
        self.sliderphiB1.setTracking(True)
        self.sliderphiB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuephiB1 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderLeftphiB1, self.sliderphiB1, self.sliderRightphiB1, self.sliderValuephiB1]):
            boxB1.addWidget(w, 2, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)

        sliderLabel = QLabel('width [km/s]')
        self.slidervthB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftvthB1 = ExtendedQLabel('0.10')
        self.sliderRightvthB1 = ExtendedQLabel('20')
        self.slidervthB1.setMinimum(0.1)
        self.slidervthB1.setMaximum(20.0)
        self.slidervthB1.setValue(20)
        self.slidervthB1.setTracking(True)
        self.slidervthB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuevthB1 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderLeftvthB1, self.slidervthB1, self.sliderRightvthB1, self.sliderValuevthB1]):
            boxB1.addWidget(w, 3, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)

        sliderLabel = QLabel('v [km/s]')
        self.slidervB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftvB1 = ExtendedQLabel('-10')
        self.sliderRightvB1 = ExtendedQLabel('10')
        self.slidervB1.setMinimum(-10)
        self.slidervB1.setMaximum(10)
        self.slidervB1.setValue(20)
        self.slidervB1.setTracking(True)
        self.slidervB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuevB1 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderLeftvB1, self.slidervB1, self.sliderRightvB1, self.sliderValuevB1]):
            boxB1.addWidget(w, 4, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)

        sliderLabel = QLabel('tau')
        self.slidertau1 = FloatSlider(Qt.Horizontal)
        self.sliderLefttau1 = ExtendedQLabel('0.0')
        self.sliderRighttau1 = ExtendedQLabel('10')
        self.slidertau1.setMinimum(0.0)
        self.slidertau1.setMaximum(10.0)
        self.slidertau1.setValue(20)
        self.slidertau1.setTracking(True)
        self.slidertau1.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuetau1 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderLefttau1, self.slidertau1, self.sliderRighttau1, self.sliderValuetau1]):
            boxB1.addWidget(w, 5, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)

        # Second component
        # 
        boxB2 = QGridLayout()
        sliderLabel = QLabel('B [G]')
        self.sliderB2 = FloatSlider(Qt.Horizontal)
        self.sliderLeftB2 = ExtendedQLabel('0')
        self.sliderRightB2 = ExtendedQLabel('1000')
        self.sliderB2.setMinimum(0.0)
        self.sliderB2.setMaximum(1000.0)
        self.sliderB2.setValue(20)
        self.sliderB2.setTracking(True)
        self.sliderB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderValueB2 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderB2, self.sliderLeftB2, self.sliderRightB2, self.sliderValueB2]):
            boxB2.addWidget(w, 0, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderB2, SIGNAL('sliderReleased()'), self.onSliderB2)
        self.connect(self.sliderLeftB2, SIGNAL('clicked()'), self.onSliderLeftB2)
        self.connect(self.sliderRightB2, SIGNAL('clicked()'), self.onSliderRightB2)

        sliderLabel = QLabel('thB [deg]')
        self.sliderthB2 = FloatSlider(Qt.Horizontal)
        self.sliderthB2.setMinimum(0)
        self.sliderthB2.setMaximum(180)
        self.sliderthB2.setValue(20)
        self.sliderthB2.setTracking(True)
        self.sliderthB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuethB2 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderthB2, self.sliderValuethB2]):
            boxB2.addWidget(w, 1, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)

        sliderLabel = QLabel('phiB [deg]')
        self.sliderphiB2 = FloatSlider(Qt.Horizontal)
        self.sliderphiB2.setMinimum(-180)
        self.sliderphiB2.setMaximum(180)
        self.sliderphiB2.setValue(20)
        self.sliderphiB2.setTracking(True)
        self.sliderphiB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuephiB2 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderphiB2, self.sliderValuephiB2]):
            boxB2.addWidget(w, 2, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)

        sliderLabel = QLabel('width [km/s]')
        self.slidervthB2 = FloatSlider(Qt.Horizontal)
        self.slidervthB2.setMinimum(0.1)
        self.slidervthB2.setMaximum(20)
        self.slidervthB2.setValue(20)
        self.slidervthB2.setTracking(True)
        self.slidervthB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuevthB2 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.slidervthB2, self.sliderValuevthB2]):
            boxB2.addWidget(w, 3, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)

        sliderLabel = QLabel('v [km/s]')
        self.slidervB2 = FloatSlider(Qt.Horizontal)
        self.slidervB2.setMinimum(-10)
        self.slidervB2.setMaximum(10)
        self.slidervB2.setValue(20)
        self.slidervB2.setTracking(True)
        self.slidervB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuevB2 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.slidervB2, self.sliderValuevB2]):
            boxB2.addWidget(w, 4, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)

        sliderLabel = QLabel('tau')
        self.slidertau2 = FloatSlider(Qt.Horizontal)
        self.slidertau2.setMinimum(0.0)
        self.slidertau2.setMaximum(10.0)
        self.slidertau2.setValue(20)
        self.slidertau2.setTracking(True)
        self.slidertau2.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuetau2 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.slidertau2, self.sliderValuetau2]):
            boxB2.addWidget(w, 5, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)


        comp1Group = QGroupBox("Component 1")
        comp1Group.setLayout(boxB1)
        vboxB1 = QVBoxLayout()
        vboxB1.addWidget(comp1Group)

        comp2Group = QGroupBox("Component 2")
        comp2Group.setLayout(boxB2)
        vboxB2 = QVBoxLayout()
        vboxB2.addWidget(comp2Group)


        vboxB = QVBoxLayout()
        vboxB.addLayout(vboxB1)
        vboxB.addLayout(vboxB2)

        # Radiative transfer
        # 
        radTranGroup = QGroupBox("Radiative transfer")
        radio1 = QRadioButton("Optically thin")
        radio2 = QRadioButton("Exact")
        radio2.setChecked(True)
        
        vboxRadTran = QHBoxLayout()
        vboxRadTran.addWidget(radio1)
        vboxRadTran.addWidget(radio2)
        vboxRadTran.addStretch(1)
        radTranGroup.setLayout(vboxRadTran)

        boundary = QGridLayout()
        self.checkAllen = QCheckBox("Use Allen")
        self.checkAllen.toggle()
        self.connect(self.checkAllen, SIGNAL('stateChanged(int)'), self.onCheckAllen)

        boundary.addWidget(self.checkAllen, 0, 0)

        title = QLabel('I0:')
        self.I0 = QLineEdit("0")
        boundary.addWidget(title, 0, 1)
        boundary.addWidget(self.I0, 0, 2)

        title = QLabel('Q0:')
        self.Q0 = QLineEdit("0")        
        boundary.addWidget(title, 0, 3)
        boundary.addWidget(self.Q0, 0, 4)

        title = QLabel('U0:')
        self.U0 = QLineEdit("0")
        boundary.addWidget(title, 0, 5)
        boundary.addWidget(self.U0, 0, 6)

        title = QLabel('V0:')
        self.V0 = QLineEdit("0")
        boundary.addWidget(title, 0, 7)
        boundary.addWidget(self.V0, 0, 8)

        title = QLabel('N slabs:')
        self.nslabs = QComboBox()
        self.nslabs.addItem("One")
        self.nslabs.addItem("Two vertical")
        self.nslabs.addItem("Two horizontal")
        self.connect(self.nslabs, SIGNAL(''), self.onChangeSlabs)
        boundary.addWidget(title, 0, 9)
        boundary.addWidget(self.nslabs, 0, 10)
        

        wave = QGridLayout()
        title = QLabel('Allen:')
        self.allen = QLineEdit("0")
        wave.addWidget(title, 0, 0)
        wave.addWidget(self.allen, 0, 1)

        title = QLabel('wl:')
        self.leftWave = QLineEdit("0")
        wave.addWidget(title, 0, 2)
        wave.addWidget(self.leftWave, 0, 3)

        title = QLabel('wr:')
        self.rightWave = QLineEdit("0")
        wave.addWidget(title, 0, 4)
        wave.addWidget(self.rightWave, 0, 5)

        title = QLabel('step:')
        self.stepWave = QLineEdit("0")
        wave.addWidget(title, 0, 6)
        wave.addWidget(self.stepWave, 0, 7)

        vboxR = QVBoxLayout()        
        vboxR.addWidget(radTranGroup)
        vboxR.addLayout(boundary)
        vboxR.addLayout(wave)


        # Radiative transfer
        # 
        boxO1 = QGridLayout()
        sliderLabel = QLabel('theta [deg]')
        self.slidertheta = FloatSlider(Qt.Horizontal)
        self.slidertheta.setMinimum(0.0)
        self.slidertheta.setMaximum(180)
        self.slidertheta.setValue(0)
        self.slidertheta.setTracking(True)
        self.slidertheta.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuertheta = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.slidertheta, self.sliderValuertheta]):
            boxO1.addWidget(w, 0, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.slidertheta, SIGNAL('sliderReleased()'), self.onSliderTheta)
        self.connect(self.sliderValuertheta, SIGNAL('clicked()'), self.onSliderValueTheta)


        sliderLabel = QLabel('phi [deg]')
        self.sliderphi = FloatSlider(Qt.Horizontal)
        self.sliderphi.setMinimum(0.1)
        self.sliderphi.setMaximum(20)
        self.sliderphi.setValue(20)
        self.sliderphi.setTracking(True)
        self.sliderphi.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuerphi = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderphi, self.sliderValuerphi]):
            boxO1.addWidget(w, 1, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderphi, SIGNAL('sliderReleased()'), self.onSliderPhi)
        self.connect(self.sliderValuerphi, SIGNAL('clicked()'), self.onSliderValuePhi)

        sliderLabel = QLabel('gamma [deg]')
        self.slidergamma = FloatSlider(Qt.Horizontal)
        self.slidergamma.setMinimum(0)
        self.slidergamma.setMaximum(180)
        self.slidergamma.setValue(0)
        self.slidergamma.setTracking(True)
        self.slidergamma.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuergamma = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.slidergamma, self.sliderValuergamma]):
            boxO1.addWidget(w, 2, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.slidergamma, SIGNAL('sliderReleased()'), self.onSliderGamma)
        self.connect(self.sliderValuergamma, SIGNAL('clicked()'), self.onSliderValueGamma)

        sliderLabel = QLabel('height [deg]')
        self.sliderheight = FloatSlider(Qt.Horizontal)
        self.sliderheight.setMinimum(0.1)
        self.sliderheight.setMaximum(20.0)
        self.sliderheight.setValue(20)
        self.sliderheight.setTracking(True)
        self.sliderheight.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuerheight = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderheight, self.sliderValuerheight]):
            boxO1.addWidget(w, 3, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderheight, SIGNAL('sliderReleased()'), self.onSliderHeight)
        self.connect(self.sliderValuerheight, SIGNAL('clicked()'), self.onSliderValueHeight)

        sliderLabel = QLabel('a')
        self.sliderDamp = FloatSlider(Qt.Horizontal)
        self.sliderDamp.setMinimum(0.0)
        self.sliderDamp.setMaximum(2.0)
        self.sliderDamp.setValue(0.0)
        self.sliderDamp.setTracking(True)
        self.sliderDamp.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuerDamp = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderDamp, self.sliderValuerDamp]):
            boxO1.addWidget(w, 4, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderDamp, SIGNAL('sliderReleased()'), self.onSliderDamp)
        self.connect(self.sliderValuerDamp, SIGNAL('clicked()'), self.onSliderValueDamp)

        sliderLabel = QLabel('ff')
        self.sliderff = FloatSlider(Qt.Horizontal)
        self.sliderff.setMinimum(0.0)
        self.sliderff.setMaximum(1.0)
        self.sliderff.setValue(0.0)
        self.sliderff.setTracking(True)
        self.sliderff.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuerff = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderff, self.sliderValuerff]):
            boxO1.addWidget(w, 5, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderff, SIGNAL('sliderReleased()'), self.onSliderff)
        self.connect(self.sliderValuerff, SIGNAL('clicked()'), self.onSliderValueff)

        sliderLabel = QLabel('S2/S1')
        self.sliderS2S1 = FloatSlider(Qt.Horizontal)
        self.sliderS2S1.setMinimum(0.0)
        self.sliderS2S1.setMaximum(2.0)
        self.sliderS2S1.setValue(0.0)
        self.sliderS2S1.setTracking(True)
        self.sliderS2S1.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuerS2S1 = ExtendedQLabel('0')
        for l, w in enumerate([sliderLabel, self.sliderS2S1, self.sliderValuerS2S1]):
            boxO1.addWidget(w, 6, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.connect(self.sliderS2S1, SIGNAL('sliderReleased()'), self.onSliderS2S1)
        self.connect(self.sliderValuerS2S1, SIGNAL('clicked()'), self.onSliderValueS2S1)

        sliderLabel = QLabel('multiplet')
        self.slidermultiplet = QSlider(Qt.Horizontal)
        self.slidermultiplet.setRange(1, 4)
        self.slidermultiplet.setValue(1)
        self.slidermultiplet.setTracking(True)
        self.slidermultiplet.setTickPosition(QSlider.TicksBothSides)
        self.sliderValuermultiplet = ExtendedQLabel('10830') 
        for l, w in enumerate([sliderLabel, self.slidermultiplet, self.sliderValuermultiplet]):
            boxO1.addWidget(w, 7, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)


        obsGroup = QGroupBox("Observation")
        obsGroup.setLayout(boxO1)
        vboxO1 = QVBoxLayout()
        vboxO1.addWidget(obsGroup)
        # vboxO1.addWidget(self.allenLabel)

        
        self.calculateButton = QPushButton("Calculate")
        self.loadButton = QPushButton("Load observation")
        self.resetButton = QPushButton("Reset observation")

        boxO2 = QVBoxLayout()
        boxO2.addWidget(self.calculateButton)
        boxO2.addWidget(self.loadButton)
        boxO2.addWidget(self.resetButton)

        calcGroup = QGroupBox("Calculate")
        calcGroup.setLayout(boxO2)
        vboxO2 = QVBoxLayout()
        vboxO2.addWidget(calcGroup)      

        vboxO = QVBoxLayout()
        vboxO.addLayout(vboxO1)
        vboxO.addLayout(vboxO2)
        
        # Final layout
        #
        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.canvas)
        hbox1.addLayout(vboxB)
        hbox1.addLayout(vboxO)

        hbox2 = QHBoxLayout()
        hbox2.addLayout(vboxR)        

        vbox = QVBoxLayout()
        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)


        # self.textbox = QLineEdit()
        # self.textbox.setMinimumWidth(200)
        # self.connect(self.textbox, SIGNAL('editingFinished ()'), self.on_draw)
        
        # self.draw_button = QPushButton("&Draw")
        # self.connect(self.draw_button, SIGNAL('clicked()'), self.on_draw)
        
        # self.grid_cb = QCheckBox("Show &Grid")
        # self.grid_cb.setChecked(False)
        # self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.on_draw)
        
        # slider_label = QLabel('Bar width (%):')
        # self.slider = QSlider(Qt.Horizontal)
        # self.slider.setRange(1, 100)
        # self.slider.setValue(20)
        # self.slider.setTracking(True)
        # self.slider.setTickPosition(QSlider.TicksBothSides)
        # self.connect(self.slider, SIGNAL('valueChanged(int)'), self.on_draw)
        
        # #
        # # Layout with box sizers
        # # 
        # hbox = QHBoxLayout()
        
        # for w in [  self.textbox, self.draw_button, self.grid_cb,
        #             slider_label, self.slider]:
        #     hbox.addWidget(w)
        #     hbox.setAlignment(w, Qt.AlignVCenter)
        
        # vbox = QVBoxLayout()
        # vbox.addWidget(self.canvas)
        # vbox.addWidget(self.mpl_toolbar)
        # vbox.addLayout(hbox)
        # vbox.addLayout(vboxB1)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("OK")
        self.statusBar().addWidget(self.status_text, 1)

    def setMultitermHe(self):
        pass

    def setMultitermNa(self):
        pass
        
    def create_menu(self):        
        """Create main menu
        
        Returns:
            TYPE: None
        """
        self.file_menu = self.menuBar().addMenu("&Multiterm")

        he_action = self.create_action("&He", slot=self.setMultitermHe, 
            shortcut="Ctrl+H", tip="Change to He I atom model")

        na_action = self.create_action("&He", slot=self.setMultitermNa, 
            shortcut="Ctrl+N", tip="Change to Na I atom model")
        
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, (he_action, na_action, None, quit_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About hazel')
        
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()