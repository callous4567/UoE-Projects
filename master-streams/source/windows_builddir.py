import windows_directories
from windows_directories import imgdir, datadir, duplimontedir
import ascii_info
import os

# Generic tools to build directories
class dirbuilder(object):
    def __init__(self):
        self.datimgdir()
        self.duplimonte()
        self.duplimonte_kmeanshtmldir()
        self.clustdir()

    # Build directory for building image/data directories
    def datimgdir(self):
        try:
            os.mkdir(imgdir)
        except:
            pass
        try:
            os.mkdir(datadir)
        except:
            pass

    # Build directory for storing duplimonte clusterings
    def duplimonte(self):
        try:
            os.mkdir(duplimontedir)
        except:
            pass
        for group in ascii_info.all_groups:
            try:
                os.mkdir(duplimontedir + "\\" + group)
            except:
                pass
        # and the full one
        try:
            os.mkdir(duplimontedir + "\\" + ascii_info.fullgroup)
        except:
            pass

    # Build directories for clustering all the duplimonte (saving the images.)
    def duplimonte_kmeanshtmldir(self):
        for group in ascii_info.all_groups:
            try:
                try:
                    os.mkdir(imgdir + "\\kmeans_html" + "\\duplimonte_kmeanshtml\\")
                except:
                    pass
                os.mkdir(imgdir + "\\kmeans_html" + "\\duplimonte_kmeanshtml\\" + group)
            except Exception as e:
                pass
        # and the full one
        try:
            os.mkdir(imgdir + "\\kmeans_html" + "\\duplimonte_kmeanshtml\\" + ascii_info.fullgroup)
        except:
            pass

    # Build cluster dir
    def clustdir(self):
        try:
            os.mkdir(windows_directories.clusterdir)
        except:
            pass
        try:
            os.mkdir(imgdir + "\\clustered")
        except:
            pass

dirbuilder()