#
# spec file for package opm-core
#

Name:           opm-upscaling
Version:        2013.03
Release:        0
Summary:        Open Porous Media - upscaling library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        %{name}-%{version}.tar.gz
BuildRequires:  blas-devel gcc-c++ gcc-gfortran lapack-devel dune-common-devel
BuildRequires:  boost-devel git suitesparse-devel cmake28 doxygen bc
BuildRequires:  tinyxml-devel dune-istl-devel superlu-devel opm-core-devel
BuildRequires:  opm-porsol-devel
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
Requires:       libopm-upscaling1 = %{version}

%description
This module provides semi-implicit pressure and transport solvers using the IMPES method.

%package -n libopm-upscaling1
Summary:        Open Porous Media - upscaling library
Group:          System/Libraries

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

%description bin
This package contains the applications for opm-upscaling

%prep
%setup -q

%build
cmake28 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF
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
