#
# spec file for package opm-upscaling
#

%define tag final

Name:           opm-upscaling
Version:        2019.04
Release:        0
Summary:        Open Porous Media - upscaling library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel dune-common-devel opm-grid-devel
BuildRequires:  git suitesparse-devel doxygen bc dune-geometry-devel
BuildRequires:  tinyxml-devel dune-istl-devel opm-common-devel
BuildRequires:  opm-material-devel ecl-devel dune-grid-devel
BuildRequires:  devtoolset-6-toolchain zlib-devel
BuildRequires:  openmpi-devel trilinos-openmpi-devel ptscotch-openmpi-devel scotch-devel opm-common-openmpi-devel opm-grid-openmpi-devel opm-material-openmpi-devel
BuildRequires:  mpich-devel trilinos-mpich-devel ptscotch-mpich-devel opm-common-mpich-devel opm-grid-mpich-devel opm-material-mpich-devel
%{?el6:BuildRequires: cmake3 boost148-devel}
%{!?el6:BuildRequires: cmake boost-devel}
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
Requires:       libopm-upscaling1 = %{version}

%description
This module provides semi-implicit pressure and transport solvers using the IMPES method.

%package -n libopm-upscaling1
Summary:        Open Porous Media - upscaling library
Group:          System/Libraries

%description -n libopm-upscaling1
This module implements single-phase and steady-state upscaling methods.

%package -n libopm-upscaling1-openmpi
Summary:        Open Porous Media - upscaling library
Group:          System/Libraries

%description -n libopm-upscaling1-openmpi
This module implements single-phase and steady-state upscaling methods.

%package -n libopm-upscaling1-mpich
Summary:        Open Porous Media - upscaling library
Group:          System/Libraries

%description -n libopm-upscaling1-mpich
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

%package openmpi-devel
Summary:        Development and header files for opm-upscaling
Group:          Development/Libraries/C and C++
Requires:       %{name} = %{version}
Requires:       blas-devel
Requires:       lapack-devel
Requires:       suitesparse-devel
Requires:       libopm-upscaling1-openmpi = %{version}

%description openmpi-devel
This package contains the development and header files for opm-upscaling

%package mpich-devel
Summary:        Development and header files for opm-upscaling
Group:          Development/Libraries/C and C++
Requires:       %{name} = %{version}
Requires:       blas-devel
Requires:       lapack-devel
Requires:       suitesparse-devel
Requires:       libopm-upscaling1-mpich = %{version}

%description mpich-devel
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

%package openmpi-bin
Summary:        Applications in opm-upscaling
Group:          Scientific
Requires:       %{name} = %{version}
Requires:       libopm-upscaling1-openmpi = %{version}
Requires:       libopm-grid1-openmpi = %{version}
Requires:       libopm-common1-openmpi = %{version}

%description openmpi-bin
This package contains the applications for opm-upscaling

%package mpich-bin
Summary:        Applications in opm-upscaling
Group:          Scientific
Requires:       %{name} = %{version}
Requires:       libopm-upscaling1-mpich = %{version}
Requires:       libopm-grid1-mpich = %{version}
Requires:       libopm-common1-mpich = %{version}

%description mpich-bin
This package contains the applications for opm-upscaling

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
scl enable devtoolset-6 bash
mkdir serial
cd serial
%{?el6:cmake3} %{!?el6:cmake} -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=/usr/include/boost148} -DINSTALL_BENCHMARKS=1 -DWITH_NATIVE=OFF ..
make %{?_smp_mflags}
make test
cd ..

mkdir openmpi
cd openmpi
%{?el6:module load openmpi-x86_64}
%{?!el6:module load mpi/openmpi-x86_64}
%{?el6:cmake3} %{!?el6:cmake} -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi -DCMAKE_INSTALL_LIBDIR=lib -DUSE_RUNPATH=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=/usr/include/boost148} -DINSTALL_BENCHMARKS=1 -DWITH_NATIVE=OFF -DZOLTAN_ROOT=/usr/lib64/openmpi -DCMAKE_CXX_FLAGS=-I/usr/include/openmpi-x86_64/trilinos -DZOLTAN_INCLUDE_DIRS=/usr/include/openmpi-x86_64/trilinos -DPTSCOTCH_ROOT=/usr/lib64/openmpi -DPTSCOTCH_INCLUDE_DIR=/usr/include/openmpi-x86_64 ..
make %{?_smp_mflags}
#make test
cd ..

mkdir mpich
cd mpich
%{?el6:module rm openmpi-x86_64}
%{?el6:module load mpich-x86_64}
%{?!el6:module rm mpi/openmpi-x86_64}
%{?!el6:module load mpi/mpich-x86_64}
%{?el6:cmake3} %{!?el6:cmake} -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/mpich -DCMAKE_INSTALL_LIBDIR=lib -DUSE_RUNPATH=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=/usr/include/boost148} -DINSTALL_BENCHMARKS=1 -DWITH_NATIVE=OFF -DZOLTAN_ROOT=/usr/lib64/mpich -DCMAKE_CXX_FLAGS=-I/usr/include/mpich-x86_64/trilinos -DZOLTAN_INCLUDE_DIRS=/usr/include/mpich-x86_64/trilinos -DPTSCOTCH_ROOT=/usr/lib64/mpich -DPTSCOTCH_INCLUDE_DIR=/usr/include/mpich-x86_64 ..
make %{?_smp_mflags}
#make test

%install
cd serial
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}
cd ..

cd openmpi
make install DESTDIR=${RPM_BUILD_ROOT}
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
cd ..

cd mpich
make install DESTDIR=${RPM_BUILD_ROOT}
mv ${RPM_BUILD_ROOT}/usr/lib64/mpich/include/* ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/
cd ..

%clean
rm -rf %{buildroot}

%post -n libopm-upscaling1 -p /sbin/ldconfig
%post -n libopm-upscaling1-openmpi -p /sbin/ldconfig
%post -n libopm-upscaling1-mpich -p /sbin/ldconfig

%postun -n libopm-upscaling1 -p /sbin/ldconfig
%postun -n libopm-upscaling1-openmpi -p /sbin/ldconfig
%postun -n libopm-upscaling1-mpich -p /sbin/ldconfig

%files
%doc README COPYING

%files doc
%{_docdir}/*

%files -n libopm-upscaling1
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files -n libopm-upscaling1-openmpi
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so.*

%files -n libopm-upscaling1-mpich
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
/usr/lib/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*
%{_datadir}/opm/cmake/Modules/*

%files openmpi-devel
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so
%{_libdir}/openmpi/lib/dunecontrol/*
%{_libdir}/openmpi/lib/pkgconfig/*
%{_includedir}/openmpi-x86_64/*
%{_libdir}/openmpi/share/cmake/*
%{_libdir}/openmpi/share/opm/cmake/Modules/*

%files mpich-devel
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so
%{_libdir}/mpich/lib/dunecontrol/*
%{_libdir}/mpich/lib/pkgconfig/*
%{_includedir}/mpich-x86_64/*
%{_libdir}/mpich/share/cmake/*
%{_libdir}/mpich/share/opm/cmake/Modules/*

%files bin
%{_bindir}/*

%files openmpi-bin
%{_libdir}/openmpi/bin/*

%files mpich-bin
%{_libdir}/mpich/bin/*
