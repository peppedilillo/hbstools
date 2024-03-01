py.install_sources(
  [
    '__init__.py',
    'bft.py',
    'bft_cwrap.py',
    'poissonfocus.py',
    'poissonfocusdes.py',
    'poissonfocusses_cwrap.py',
  ],
  subdir: 'hbstools' / 'triggers'
)

incdir = include_directories('c' / 'include')
sources = ['c' / 'poissonfocus.c', 'c' / 'poissonfocusses.c', 'c' / 'bft.c']
shared_library('lib-pfocus',
    sources,
    include_directories : incdir,
    install: true,
    name_prefix: '',
    install_dir: py.get_install_dir() / 'hbstools' / 'triggers' / '.sharedlibs'
)