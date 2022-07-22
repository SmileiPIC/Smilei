import sys
import os
import getopt
try:
    import sysconfig
except:
    from distutils import sysconfig

valid_opts = ['prefix', 'exec-prefix', 'includes', 'libs', 'cflags',
              'ldflags', 'help']

def exit_with_usage(code=1):
    try:
        print("Usage: %s [%s]" % (sys.argv[0],
                                                '|'.join('--'+opt for opt in valid_opts)))
    except:
        print("Usage: {0} [{1}]".format(
            sys.argv[0], '|'.join('--'+opt for opt in valid_opts)))
    sys.exit(code)

try:
    opts, args = getopt.getopt(sys.argv[1:], '', valid_opts)
except getopt.error:
    exit_with_usage()

if not opts:
    exit_with_usage()

pyver = sysconfig.get_config_var('VERSION')
getvar = sysconfig.get_config_var

opt_flags = [flag for (flag, val) in opts]

if '--help' in opt_flags:
    exit_with_usage(code=0)

for opt in opt_flags:
    if opt == '--prefix':
        try:
            print(sysconfig.PREFIX)
        except:
            print(sysconfig.get_config_var('prefix'))
    
    elif opt == '--exec-prefix':
        try:
            print(sysconfig.EXEC_PREFIX)
        except:
            print(sysconfig.get_config_var('exec_prefix'))
    
    elif opt in ('--includes', '--cflags'):
        try:
            flags = ['-I' + sysconfig.get_python_inc(),
                     '-I' + sysconfig.get_python_inc(plat_specific=True)]
        except:
            flags = ['-I' + sysconfig.get_path('include'),
                     '-I' + sysconfig.get_path('platinclude')]
        
        try:
            import numpy
            flags += ['-I' + numpy.get_include() + ' -DSMILEI_USE_NUMPY -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']
        except:
            pass
        
        if opt == '--cflags':
            flags.extend(getvar('CFLAGS').split())
        print(' '.join(flags))

    elif opt in ('--libs', '--ldflags'):
        try:
            libs = ['-lpython' + pyver + sys.abiflags]
        except:
            libs = ['-lpython' + pyver]
        libs += getvar('LIBS').split()
        libs += getvar('SYSLIBS').split()
        # add the prefix/lib/pythonX.Y/config dir, but only if there is no
        # shared library in prefix/lib/.
        if opt == '--ldflags':
            libs.insert(0, '-L' + getvar('LIBDIR'))
            if not getvar('Py_ENABLE_SHARED'):
                libs.insert(0, '-L' + getvar('LIBPL'))
            if not getvar('PYTHONFRAMEWORK'):
                libs.extend(getvar('LINKFORSHARED').split())
        print(' '.join(libs))

