# -*- mode: python ; coding: utf-8 -*-

import periodictable
periodictable_path = periodictable.__path__[0]

import os
cwd = os.getcwd()

added_files = [
        (os.path.join(periodictable_path,'xsf/*'), 'periodictable/xsf'),
        (os.path.join(cwd,'GUI/mainwindow.ui'),'GUI')
        ]


block_cipher = None
a = Analysis(['main.py'],
             pathex=[cwd],
             binaries=[],
             datas=added_files,
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='XFNTR Analyzer',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='XFNTR Analyzer')
app = BUNDLE(coll,
             name='XFNTR Analyzer.app',
             icon=cwd+'/images/TaiChi.icns',
             bundle_identifier=None)
