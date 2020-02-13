# Developter log
## 1. Add command script: xfntr
Use "entry_points" keyword in setup() to add command "xfntr" in ocmmand lien. There might be path problems.

```
entry_points = { # create scripts and add to sys.PATH
        'console_scripts':[
            'xfntr1 = xfntr.main:main'
        ],
        'gui_scripts': [
            'xfntr = xfntr.main:main'
        ]
    }
```

## 2. Solve path issue

command line might not be albe to find the right path where the module is. 

In main.py, added the following so it can find module mainwindow.py
```
import os
# Use absolute path instead of relative path ('./') to avoid trouble when installed by pip
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)
```

In mainwindow.py, added the following
```
import os
# Use absolute path instead of relative path ('./') to avoid trouble when installed by pip
dir_path = os.path.dirname(os.path.realpath(__file__))
```
and used the absolute file path for GUI files:
```
UI_path = dir_path + '/GUI/'
```
