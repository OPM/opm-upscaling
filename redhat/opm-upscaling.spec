#
# spec file for package opm-upscaling
#

%define tag rc2

Name:           opm-upscaling
Version:        2013.10
Release:        0
Summary:        Open Porous Media - upscaling library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel dune-common-devel
BuildRequires:  git suitesparse-devel cmake28 doxygen bc
BuildRequires:  tinyxml-devel dune-istl-devel superlu-devel opm-core-devel
BuildRequires:  opm-porsol-devel opm-material-devel ert.ecl-devel
%{?el5:BuildRequires: gcc44 gcc44-gfortran gcc44-c++}
%{!?el5:BuildRequires: gcc gcc-gfortran gcc-c++}
%{?el5:BuildRequires: boost141-devel}
%{!?el5:BuildRequires: boost-devel}
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
Requires:       libopm-upscaling1 = %{version}

%description
This module provides semi-implicit pressure and transport solvers using the IMPES method.

%package -n libopm-upscaling1
Summary:        Open Porous Media - upscaling library
Group:          System/Libraries
%{?el5:BuildArch: %{_arch}}

%description -n libopm-upscaling1
This module implements single-phase and steady-state upscaling methods.

%package devel
Summary:        Development and header files for opm-upscaling
Group:          Development/Libraries/C and C++
Requires:       %{name} = %{version}
Requires:       blas-devel
Requires:       lapack-devel
Requires:       suitesparse-devel
Requires:       libopm-upscaling1 = %{version}
%{?el5:BuildArch: %{_arch}}

%description devel
This package contains the development and header files for opm-upscaling

%package doc
Summary:        Documentation files for opm-upscaling
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-upscaling

%package bin
Summary:        Applications in opm-upscaling
Group:          Scientific
Requires:       %{name} = %{version}
Requires:       libopm-upscaling1 = %{version}
%{?el5:BuildArch: %{_arch}}

%description bin
This package contains the applications for opm-upscaling

%{?el5:
%package debuginfo
Summary:        Debug info in opm-upscaling
Group:          Scientific
Requires:       %{name} = %{version}
Requires:       libopm-upscaling1 = %{version}, opm-upscaling-bin = %{version}
BuildArch: 	%{_arch}

%description debuginfo
This package contains the debug symbols for opm-upscaling
}

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
cmake28 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF %{?el5:-DCMAKE_CXX_COMPILER=g++44 -DCMAKE_C_COMPILER=gcc44 -DCMAKE_Fortran_COMPILER=gfortran44 -DBOOST_LIBRARYDIR=%{_libdir}/boost141 -DBOOST_INCLUDEDIR=/usr/include/boost141}
make

%install
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}

%clean
rm -rf %{buildroot}

%post -n libopm-upscaling1 -p /sbin/ldconfig

%postun -n libopm-upscaling1 -p /sbin/ldconfig

%files
%doc README

%files doc
%{_docdir}/*

%files -n libopm-upscaling1
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*

%files bin
%{_bindir}/*

%{?el5:
%files debuginfo
/usr/lib/debug/%{_libdir}/*.so*.debug
}
