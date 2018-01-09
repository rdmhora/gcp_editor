from PIL.ExifTags import TAGS, GPSTAGS
from PIL import Image
import numpy as np
import utm  # library used to convert coordinate system
import os
os.system("pyuic5 -x gcp_editor_gui.ui -o gcp_editor_gui.py")
os.system("pyuic5 -x saved.ui -o saved.py")
os.system("pyuic5 -x ops.ui -o ops.py")

from os.path import basename
from PyQt5 import uic, QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog, QTableWidgetItem, QTableView, QGraphicsItem, QGraphicsView
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from gcp_editor_gui import Ui_GCPEditor
from saved import Ui_Saved
from ops import Ui_Ops


import sys

gcp_file = []
imageList = 0
gcpTable = 0
images = []
zoom = 0

class CheckExistence:

    def exif_check(self, imagess):
        img = []
        for i in range(0, (len(imagess) - 1)):
            image = Image.open(imagess[i])
            img.append(image)
            exif = img[i]._getexif()
            ret = {}
            gps = {}
            for tag, value in exif.items():
                decoded = TAGS.get(tag, tag)
                ret[decoded] = value
                # Getting GPS lat, long
                if decoded == "GPSInfo":
                    gps_data = {}
                    for t in value:
                        sub_decoded = GPSTAGS.get(t, t)
                        gps_data[sub_decoded] = value[t]

                    gps[decoded] = gps_data
            test = {}
            if gps == test:
                return 0
            else:
                return 1

class GetIfExist:
    def _get_if_exist(self, data, key):
        if key in data:
            return data[key]

        return None

# Get N and E from Exif info

class GetNE:
    def get_N_E(self, gps_data):
        """Returns the latitude and longitude, if available, from the provided exif_data (obtained through get_exif_data above)"""
        lat = None
        lon = None
        get = GetIfExist()
        conv = ConvertToDegrees()
        if "GPSInfo" in gps_data:
            gps_info = gps_data["GPSInfo"]

            gps_latitude = get._get_if_exist(gps_info, "GPSLatitude")
            gps_latitude_ref = get._get_if_exist(gps_info, 'GPSLatitudeRef')
            gps_longitude = get._get_if_exist(gps_info, 'GPSLongitude')
            gps_longitude_ref = get._get_if_exist(gps_info, 'GPSLongitudeRef')

            if gps_latitude and gps_latitude_ref and gps_longitude and gps_longitude_ref:
                lat = conv._convert_to_degress(gps_latitude)
                if gps_latitude_ref != "N":
                    lat = 0 - lat

                lon = conv._convert_to_degress(gps_longitude)
                if gps_longitude_ref != "E":
                    lon = 0 - lon
        E, N, fuse, hem = utm.from_latlon(lat, lon)
        return N, E, fuse


# Read Gcp File
class ReadGcp:

    def read_gcp(self, x):
        file = open(x)
        a = file.readlines()
        b = []
        for i in range(len(a)):
            b.append(a[i].split())
        return b

# Convert lat long to UTM
class ConvertToDegrees:

    def _convert_to_degress(self, value):
        """Helper function to convert the GPS coordinates stored in the EXIF to degress in float format"""
        d0 = value[0][0]
        d1 = value[0][1]
        d = float(d0) / float(d1)

        m0 = value[1][0]
        m1 = value[1][1]
        m = float(m0) / float(m1)

        s0 = value[2][0]
        s1 = value[2][1]
        s = float(s0) / float(s1)

        return d + (m / 60.0) + (s / 3600.0)

# Get focal Length from Exif File

class FocalLen:
    def focal_leng(self, ret):
        # Focal Length in meters
        foc_leng = ret['FocalLength']
        foc_leng = 0.01 * foc_leng[0] / foc_leng[1]
        return foc_leng

# Calculate GSD from exif file

class Gsd:
    def gsd(self, ret, alt, fl, e, n):
        # Pixel size in meters

        pix_siz = ret['XResolution']

        gsde = float(e) * float(fl) / (float(pix_siz[0]) * float(alt))
        gsdn = float(n) * float(fl) / (float(pix_siz[0]) * float(alt))
        return gsde, gsdn

# Open Exif information

class OpenExif:
    def open_exif(self, img):

        exif = img._getexif()
        ret = {}
        gps = {}
        for tag, value in exif.items():
            decoded = TAGS.get(tag, tag)
            ret[decoded] = value
            # Getting GPS lat, long
            if decoded == "GPSInfo":
                gps_data = {}
                for t in value:
                    sub_decoded = GPSTAGS.get(t, t)
                    gps_data[sub_decoded] = value[t]

                gps[decoded] = gps_data
        return ret, gps

# Read XMP info to get flight information (pitch, row, yaw, relative altitue) better for DJI images

class Xmp:
    def xmp_info(self, f):
        # Getting Altitude and Flight Yaw from XMP metadata
        with open(f, "rb") as fin:
            img = fin.read()
        imgAsString = str(img)

        rel_alt = imgAsString.find('RelativeAltitude="')
        fli_yaw = imgAsString.find('FlightYawDegree="')
        alt = ' '
        yaw = ' '
        if rel_alt != -1:
            alt = imgAsString[(rel_alt + 18):(rel_alt + 24)]
            alt = float(alt)
        else:
            alt = 'no_xmp'
        if fli_yaw != -1:
            yaw = imgAsString[(fli_yaw + 17):(fli_yaw + 23)]
            yaw = float(yaw)
        else:
            yaw = 'no_xmp'

        return alt, yaw

# Calculate GSD, N, E, yaw for each image

class GetInfo:

    def info_extrat(self, file):
        img = []
        e = []
        n = []
        ret = []
        gps = []
        alt = []
        yaw = []
        gs = []
        Ngps = []
        Egps = []
        exif = OpenExif()
        xmp = Xmp()
        gsd = Gsd()
        fl = FocalLen()
        ne = GetNE()
        for i in range(0, (len(file) - 1)):
            img.append(file[i])
            image = Image.open(img[i])
            em, nm = np.size(image)
            e.append(em)
            n.append(nm)
            retm, gpsm = exif.open_exif(image)
            ret.append(retm)
            gps.append(gpsm)
            altm, yawm = xmp.xmp_info(img[i])
            alt.append(altm)
            yaw.append(yawm)
            gsm = gsd.gsd(ret[i], alt[i], fl.focal_leng(ret[i]), e[i], n[i])
            gs.append(gsm)
            Ngpsm, Egpsm, fuse = ne.get_N_E(gps[i])
            Ngps.append(Ngpsm)
            Egps.append(Egpsm)

        return gs, Egps, Ngps, yaw

# Calculate image coordinate for each GCP that is present on image
class ImageCoord:

    def coord_distribute(self, f, gsd, Eexif, Nexif, E, N, yaw):

        img = Image.open(f)
        (e, n) = np.size(img)
        new = []
        img_coor = []
        test = 1

        if abs(yaw) < 90:
            Eexif = (Eexif * round(np.cos(np.pi * yaw / 180), 2) - Nexif * round(np.sin(np.pi * yaw / 180), 2))
            Nexif = (Eexif * round(np.sin(np.pi * yaw / 180), 2) + Nexif * round(np.cos(np.pi * yaw / 180), 2))
            E = (E * round(np.cos(np.pi * yaw / 180), 2) - N * round(np.sin(np.pi * yaw / 180), 2))
            N = (E * round(np.sin(np.pi * yaw / 180), 2) + N * round(np.cos(np.pi * yaw / 180), 2))
            new.append([Eexif - gsd[0] * e / 2, Nexif - gsd[0] * n / 2])
            new.append([Eexif + gsd[0] * e / 2, Nexif + gsd[0] * n / 2])

            if new[0][0] < E < new[1][0] and new[0][1] < N < new[1][1]:
                img_coor.append([(e / 2 - (Eexif - E) / round(gsd[0], 3)), (n / 2 + (Nexif - N) / round(gsd[1], 3))])
            else:
                test = 0

        elif abs(yaw) > 90:
            yaw = 180 - abs(yaw)
            Eexif = (Eexif * round(np.cos(np.pi * yaw / 180), 2) - Nexif * round(np.sin(np.pi * yaw / 180), 2))
            Nexif = (Eexif * round(np.sin(np.pi * yaw / 180), 2) + Nexif * round(np.cos(np.pi * yaw / 180), 2))
            E = (E * round(np.cos(np.pi * yaw / 180), 2) - N * round(np.sin(np.pi * yaw / 180), 2))
            N = (E * round(np.sin(np.pi * yaw / 180), 2) + N * round(np.cos(np.pi * yaw / 180), 2))
            new.append([Eexif - gsd[0] * e / 2, Nexif - gsd[0] * n / 2])
            new.append([Eexif + gsd[0] * e / 2, Nexif + gsd[0] * n / 2])

            if new[0][0] < E < new[1][0] and new[0][1] < N < new[1][1]:
                img_coor.append([(e / 2 + (Eexif - E) / (round(gsd[0], 2))), (n / 2 - (Nexif - N) / (round(gsd[1], 2)))])
            else:
                test = 0

        return img_coor, test


class MyScene(QtWidgets.QGraphicsScene):

    def mousePressEvent(self, e):
        from Gcp_Editor_v1 import GCPEditor
        global gcp_file
        global imageList
        global gcpTable
        global images

        #try:
        for i in range(1,len(gcp_file)-1):
                if imageList.item(imageList.currentRow()).text() == gcp_file[i][5] and \
                                    gcpTable.item(gcpTable.currentRow(), 0).text() == gcp_file[i][0]:

                    gcp_file[i][4][0][0] = e.scenePos().x()
                    gcp_file[i][4][0][1] = e.scenePos().y()
                    pen = QPen(Qt.blue)
                    pen.setCosmetic(True)
                    self.addLine(e.scenePos().x(), e.scenePos().y() - 20, e.scenePos().x(), e.scenePos().y() + 20, pen)
                    self.addLine(e.scenePos().x()-20, e.scenePos().y(), e.scenePos().x() + 20, e.scenePos().y(), pen)
        self.update()

        #except:
         #   print('nothing')



class MyView(QtWidgets.QGraphicsView):

    def wheelEvent(self, event):
        try:
            global images
            global zoom
            self._zoom = 0
            img = Image.open(images[0])
            (n, e) = img.size
            if event.angleDelta().y() / 2800 > 0:
                factor = 1.25
                self._zoom += 1
                zoom = self._zoom + zoom
            else:
                factor = 0.8
                self._zoom -= 1
                zoom = self._zoom + zoom
            if self._zoom >= 0:
                self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
                self.scale(factor, factor)
            elif self._zoom < 0:
                self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
                self.scale(factor, factor)

            if zoom < 0:
                self.fitInView(QRectF(0, 0, n, e), Qt.KeepAspectRatio)
                zoom = 0
        except:
            print('nothing')


class GCPEditor(QtWidgets.QMainWindow):

    def __init__(self):
        super(GCPEditor, self).__init__()
        self.ui = Ui_GCPEditor()
        self.ui.setupUi(self)
        self.ui.gcpTable.setColumnCount(4)
        self.ui.gcpTable.setRowCount(1)
        self.ui.gcpTable.setHorizontalHeaderLabels(["Name", "E", "N", "Z"])
        self.items = self.ui.gcpTable.setHorizontalHeaderLabels(["Name", "E", "N", "Z"])
        self.ui.gcpTable.setSelectionBehavior(QTableView.SelectRows)
        self.ui.gcpTable.clicked.connect(self.gcpListClicked)
        self.ui.imageList.clicked.connect(self.imageListClicked)
        self.ui.saveButton.clicked.connect(self.save)
        self.ui.addImagesButton.clicked.connect(self.addImages)
        self.ui.addGcpFile.clicked.connect(self.addGcp)
        self.ui.pushButton.clicked.connect(self.run)
        self.ui.imageScene = MyScene(self.ui.centralwidget)
        self.ui.imageScene.setObjectName("imageScene")
        self.ui.imageView = MyView(self.ui.centralwidget)
        self.ui.imageView.setGeometry(QtCore.QRect(20, 320, 671, 421))
        self.ui.imageView.setFrameShape(QtWidgets.QGraphicsView.StyledPanel)
        self.ui.imageView.setFrameShadow(QtWidgets.QGraphicsView.Raised)
        self.ui.imageView.setObjectName("imageView")
        self.ui.imageView.setScene(self.ui.imageScene)

        global imageList
        global gcpTable

        imageList = self.ui.imageList
        gcpTable = self.ui.gcpTable

        self.images = 0

        self.alt = 0
        self.y = 0
        self.zoom = 1

    def addImages(self):
        global images
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.ExistingFiles)
        try:
            if dlg.exec_():

                self.images = dlg.selectedFiles()
            check = CheckExistence()
            exif_test = check.exif_check(self.images)
            # If aexist Exif list the gcp coordinates at table

            if exif_test is 1:
                self.ui.imageLabel.setText(str(len(self.images)) + ' Images Added')
                images = self.images
            else:
                self.images = 0

                self.addImg = QtWidgets.QMainWindow()
                self.addUi = Ui_Ops()
                self.addUi.setupUi(self.addImg)
                self.addImg.show()
        except:

            print('')


    def addGcp(self):
        global gcp_file
        dlg = QFileDialog()
        try:
            if dlg.exec_():
                gcp_file = dlg.selectedFiles()
            gcp = open(gcp_file[0])
            self.ui.gcpLabel.setText(str(len(gcp.readlines())) + ' GCPs Added')
        except:
            print('nothing')

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def rowCount(self, parent):
        return 1

    def columnCount(self, parent):
        return len(self.items)

    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        elif role != Qt.DisplayRole:
            return QVariant()

        column = index.column()
        if column < len(self.items):
            return QVariant(self.items[column])
        else:
            return QVariant()

            #######################################################################################################

    def imageListClicked(self, clickedIndex):
        global gcp_file
        index = clickedIndex.row()
        self.ui.imageScene.clear()
        for i in range(0, len(self.images) - 1):
                if self.ui.imageList.item(index).text() == basename(self.images[i]):
                    for j in range(1, len(gcp_file) - 1):
                        if self.ui.gcpTable.item(self.ui.gcpTable.currentRow(), 0).text() == gcp_file[j][0] and \
                                     self.ui.imageList.item(index).text() == gcp_file[j][5]:
                            img = Image.open(self.images[i])
                            (n, e) = img.size

                            pixMap = QPixmap(self.images[i])
                            self.ui.imageScene.addPixmap(pixMap)
                            self.ui.imageView.fitInView(QRectF(0, 0, n, e), Qt.KeepAspectRatio)
                            pen = QPen(Qt.green)
                            pen.setCosmetic(True)
                            self.ui.imageScene.addLine(gcp_file[j][4][0][0], gcp_file[j][4][0][1] - 20,
                                                              gcp_file[j][4][0][0], gcp_file[j][4][0][1] + 20, pen)
                            self.ui.imageScene.addLine(gcp_file[j][4][0][0] - 20, gcp_file[j][4][0][1],
                                                              gcp_file[j][4][0][0] + 20, gcp_file[j][4][0][1], pen)
                            self.ui.imageScene.update()

    def gcpListClicked(self, clickedIndex):
        global gcp_file

        row = clickedIndex.row()
        model = clickedIndex.model()
        column = clickedIndex.column()
        self.ui.imageList.clear()
        # add at image List the images that contain gcp
        for i in range(1, (len(gcp_file) - 1)):
            if gcp_file[i][0] == self.ui.gcpTable.item(row, 0).text():
                self.ui.imageList.addItems([gcp_file[i][5]])
        self.ui.imageList.setCurrentRow(0)

    def run(self):
        global gcp_file
        try:
            read = ReadGcp()
            check = CheckExistence()
            gcp = read.read_gcp(gcp_file[0])
            # Check is images have EXIF
            exif_test = check.exif_check(self.images)
            # If aexist Exif list the gcp coordinates at table

            if exif_test is 1:
                self.ui.gcpTable.setColumnCount(len(gcp[0]))
                self.ui.gcpTable.setRowCount(len(gcp))
                for i in range(0, len(gcp[0])):
                    for j in range(0, len(gcp)):
                        self.ui.gcpTable.setItem(j, i, QTableWidgetItem(gcp[j][i]))
                self.ui.gcpTable.selectRow(0)


                # separate images that have each GCP
                self.alt = []
                self.y = []
                xmp = Xmp()
                info = GetInfo()
                imgcoord = ImageCoord()
                for i in range(0, len(self.images) - 1):
                    aa, yy = xmp.xmp_info(self.images[i])
                    self.alt.append(aa)
                    self.y.append(yy)
                gs1, e1, n1, y1 = info.info_extrat(self.images)
                gcpname = []
                # Automatic cood finder
                for i in range(0, len(gcp)):
                    gcpname.append(gcp[i][0])
                    E = float(gcp[i][1])
                    N = float(gcp[i][2])
                    Z = float(gcp[i][3])
                    for j in range(0, (len(self.images) - 1)):
                        coord, test = imgcoord.coord_distribute(self.images[j], gs1[j], e1[j], n1[j], E, N, y1[j])
                        img = Image.open(self.images[j])
                        e, n = np.size(img)
                        if test is 1 and e > coord[0][0] > 0 and n > coord[0][1] > 0:
                            gcp_file.append([gcpname[i], E, N, Z, coord, basename(self.images[j])])

            else:
                print('')


        except:
            print('Create a dialog asking for images or cgp file')

    def save(self):
        global gcp_file
        gcp_file.remove(gcp_file[0])

        try:
            thefile = open('gcp_list.txt', 'w')
            for item in gcp_file:
                thefile.write("%s\t" % str(item[1]) + "%s\t" % str(item[2]) + "%s\t" % str(item[3]) +
                                  "%s\t" % str(int(item[4][0][0])) + "%s\t" % str(int(item[4][0][1])) + "%s\n" % str(item[5]))
            thefile.close()
            self.save = QtWidgets.QMainWindow()
            self.saveUi = Ui_Saved()
            self.saveUi.setupUi(self.save)
            self.save.show()

        except:
            print('')


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = GCPEditor()
    window.show()
    sys.exit(app.exec_())

