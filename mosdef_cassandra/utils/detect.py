import shutil


def detect_cassandra_binaries():

    cassandra_exec_names = [
        "cassandra_gfortran_openMP.exe",
        "cassandra_pgfortran_openMP.exe",
        "cassandra_intel_openMP.exe",
        "cassandra.exe",
        "cassandra_gfortran.exe",
        "cassandra_pgfortran.exe",
    ]

    py2_exec_names = ["python2", "python2.7"]

    for name in cassandra_exec_names:
        cassandra = shutil.which(name)
        if cassandra is not None:
            break

    fraglib_setup = shutil.which("library_setup.py")

    if cassandra is None or fraglib_setup is None:
        raise ValueError(
            "Error detecting cassandra. Both 'cassandra_*.exe' and "
            "'library_setup.py' must be in your PATH"
        )

    for name in py2_exec_names:
        py2 = shutil.which(name)
        if py2 is not None:
            break
    if py2 is None:
        raise ValueError(
            "Error detecting python2. library_setup.py requires python2"
        )

    print("Using the following executables for Cassandra:")
    print("Python: {}".format(py2))
    print("library_setup: {}".format(fraglib_setup))
    print("Cassandra: {}".format(cassandra))

    return py2, fraglib_setup, cassandra
