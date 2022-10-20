import ascii_info
import windows_directories
import os
import shutil

# Literally just copy the basic asciiname so that finetune has something to work with
os.chdir(windows_directories.datadir)
shutil.copyfile(os.path.join(windows_directories.datadir, ascii_info.asciiname),
                os.path.join(windows_directories.datadir, ascii_info.finetune_asciiname))