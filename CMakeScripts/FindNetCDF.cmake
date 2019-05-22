#[[ FindNetCDF.cmake

This module finds NetCDF-C and NetCDF-F. It uses nc-config and nf-config to
get HINTS for the find_xxxx's that are used to find the files/directories 
listed below.

If a file or directory cannot be found, the user should add the appropriate
directories to CMAKE_PREFIX_PATH.

Resulting variables:
    NETCDF_F_LIBRARY:       Path to libnetcdff.so
    NETCDF_C_LIBRARY:       Path to libnetcdf.so
    NETCDF_C_INCLUDE_DIR:   Path to the directory containing netcdf.h
    NETCDF_F90_INCLUDE_DIR: Path to the directory containing netcdf.mod
    NETCDF_F77_INCLUDE_DIR: Path to the directory containing netcdf.inc

    NETCDF_LIBRARIES:       Paths to all of NetCDF's libraries
    NETCDF_INCLUDE_DIRS:    Paths to all of NetCDF's include directories.

]]


# Find the nc-config and nf-config programs
find_program(NC_CONFIG NAMES "nc-config" DOC "Location of nc-config utility")
find_program(NF_CONFIG NAMES "nf-config" DOC "Location of nf-config utility")
if(NOT NC_CONFIG)
    message(WARNING "\nCouldn't find nc-config. Is it in your PATH?\n")
endif()
if(NOT NF_CONFIG)
    message(WARNING "\nCouldn't find nf-config. Is it in your PATH?\n")
endif()

# A function to call nx-config with an argument, and append the resulting path to a list
function(inspect_netcdf_config VAR NX_CONFIG ARG)
    execute_process(
        COMMAND ${NX_CONFIG} ${ARG}
        RESULT_VARIABLE NX_CONFIG_RET   
        OUTPUT_VARIABLE NX_CONFIG_OUTPUT
        ERROR_VARIABLE  NX_CONFIG_STDERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(EXISTS "${NX_CONFIG_OUTPUT}")
        list(APPEND ${VAR} ${NX_CONFIG_OUTPUT})
        set(${VAR} ${${VAR}} PARENT_SCOPE)
    endif()
endfunction()

# Determine HINTS for netcdf.h
set(NC_INC_HINTS "")
inspect_netcdf_config(NC_INC_HINTS ${NC_CONFIG} --includedir)
inspect_netcdf_config(NC_INC_HINTS ${NC_CONFIG} --prefix)
# Find netcdf.h
find_path(NETCDF_C_INCLUDE_DIR
    netcdf.h
    DOC "Directory containing \"netcdf.h\""
    HINTS ${NC_INC_HINTS}
    PATH_SUFFIXES
        "include"
)

# Determine HINTS for netcdf.mod
set(NF_INC_HINTS "")
inspect_netcdf_config(NF_INC_HINTS ${NF_CONFIG} --includedir)
inspect_netcdf_config(NF_INC_HINTS ${NF_CONFIG} --prefix)
# Find netcdf.mod
find_path(NETCDF_F90_INCLUDE_DIR
    netcdf.mod
    DOC "Directory containing \"netcdf.mod\""
    HINTS ${NF_INC_HINTS}
    PATH_SUFFIXES
        "include"
        "mod"
        "module"
)
# Find netcdf.inc
find_path(NETCDF_F77_INCLUDE_DIR
    netcdf.inc
    DOC "Directory containing \"netcdf.inc\""
    HINTS ${NF_INC_HINTS}
    PATH_SUFFIXES
        "include"
        "mod"
        "module"
)

# Determine HINTS for NetCDF-C's library
set(NC_LIBDIR_HINTS "")
inspect_netcdf_config(NC_LIBDIR_HINTS ${NC_CONFIG} --libdir)
inspect_netcdf_config(NC_LIBDIR_HINTS ${NC_CONFIG} --prefix)
# Find libnetcdf.so
find_library(NETCDF_C_LIBRARY
    netcdf
    DOC "Path to \"libnetcdf\""
    HINTS ${NC_LIBDIR_HINTS}
    PATH_SUFFIXES
        "lib"
)

# Determine HINTS for NetCDF-F's library
set(NF_LIBDIR_HINTS "")
inspect_netcdf_config(NF_LIBDIR_HINTS ${NF_CONFIG} --libdir)
inspect_netcdf_config(NF_LIBDIR_HINTS ${NF_CONFIG} --prefix)
# Find libnetcdff.so
find_library(NETCDF_F_LIBRARY
    netcdff
    DOC "Path to \"libnetcdff\""
    HINTS ${NF_LIBDIR_HINTS}
    PATH_SUFFIXES
        "lib"
)

# Make a readable error message 
set(NetCDF_ERRMSG "\nCouldn't find the following directories files:")
if(NOT NETCDF_F_LIBRARY)
    set(NetCDF_ERRMSG "${NetCDF_ERRMSG}
    NETCDF_F_LIBRARY: Path to \"libnetcdff.so\"")
endif()
if(NOT NETCDF_C_LIBRARY)
    set(NetCDF_ERRMSG "${NetCDF_ERRMSG}
    NETCDF_C_LIBRARY: Path to \"libnetcdf.so\"")
endif()
if(NOT NETCDF_C_INCLUDE_DIR)
    set(NetCDF_ERRMSG "${NetCDF_ERRMSG}
    NETCDF_C_INCLUDE_DIR: Directory containing \"netcdf.h\"")
endif()
if(NOT NETCDF_F90_INCLUDE_DIR)
    set(NetCDF_ERRMSG "${NetCDF_ERRMSG}
    NETCDF_F90_INCLUDE_DIR: Directory containing \"netcdf.mod\"")
endif()
if(NOT NETCDF_F77_INCLUDE_DIR)
    set(NetCDF_ERRMSG "${NetCDF_ERRMSG}
    NETCDF_F77_INCLUDE_DIR: Directory containing \"netcdf.inc\"")
endif()
set(NetCDF_ERRMSG "${NetCDF_ERRMSG}\n\n Set CMAKE_PREFIX_PATH to a semicolon-separated list of directories with these files.\n\n")

# Conform to the find_package standards
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
    REQUIRED_VARS 
        NETCDF_F_LIBRARY 
        NETCDF_C_LIBRARY
        NETCDF_C_INCLUDE_DIR 
        NETCDF_F90_INCLUDE_DIR 
        NETCDF_F77_INCLUDE_DIR 
    FAIL_MESSAGE "${NetCDF_ERRMSG}"
)

# Set NETCDF_LIBRARIES NETCDF_INCLUDE_DIRS
set(NETCDF_LIBRARIES ${NETCDF_F_LIBRARY} ${NETCDF_C_LIBRARY})
set(NETCDF_INCLUDE_DIRS ${NETCDF_F90_INCLUDE_DIR} ${NETCDF_F77_INCLUDE_DIR} ${NETCDF_C_INCLUDE_DIR})
