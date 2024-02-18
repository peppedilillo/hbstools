py.install_sources(
  [
    '__init__.py',
    'poissonfocus.py',
    'poissonfocusdes.py',
    'cpoissonfocusses.py',
  ],
  subdir: 'hbstools/triggers'
)

incdir = include_directories('c/include')
sources = ['c/poissonfocus.c', 'c/poissonfocusses.c', 'c/bft.c']
shared_library('c_pfocus',
    sources,
    include_directories : incdir,
    install: true,
    install_dir: py.get_install_dir() / 'hbstools/triggers/.sharedlibs'
)