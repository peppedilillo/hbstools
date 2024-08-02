project(
    'Hermes Burst Search Tools', 'c',
    default_options: [
        'buildtype=debugoptimized',
        'b_ndebug=if-release',
        'c_std=c99',
    ],
)
add_global_arguments(language : 'c')
add_project_link_arguments('-lm', language : 'c')

py = import('python').find_installation()

subdir('hbstools')
subdir('mercury')
