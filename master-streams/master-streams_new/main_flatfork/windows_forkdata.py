import ascii_info_new
import windows_directories_new
import os
import shutil

# Literally just copy the basic asciiname so that flatfork has something to work with
os.chdir(windows_directories_new.datadir)
shutil.copyfile(os.path.join(windows_directories_new.datadir, ascii_info_new.asciiname),
                os.path.join(windows_directories_new.datadir, ascii_info_new.flatfork_asciiname))