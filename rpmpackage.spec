# These xxx markers are to be replaced by git_build_rpm
Name:           xxx
Version:        xxx

Release:        1%{?dist}
Summary:        Dot (or grid) finder for camera calibrations

License:        LGPL-2.1+
URL:            https://github.com/dkogan/mrgingham/
Source0:        https://github.com/dkogan/mrgingham/archive/%{version}.tar.gz#/%{name}-%{version}.tar.gz

BuildRequires:  opencv-devel >= 3.2
BuildRequires:  boost-devel
BuildRequires:  chrpath


# to build the manpages I need to run 'mrgingham-observe-pixel-uncertainty
# --help'
BuildRequires:  python36

# for the python interface
BuildRequires: python36-numpy
BuildRequires: python36-devel
BuildRequires: python36-libs
# need this too for the /usr/include/numpy link. It points to the same place as
# the link in python36-numpy would point to, if they bothered to make one
BuildRequires: python2-numpy

Conflicts: mrgingham-tools <= 1.10

%description
Library to find a grid of points; used for calibration routines

%package devel
Requires:       %{name}%{_isa} = %{version}-%{release}
Summary:        Development files for mrgingham

%description devel
Headers and libraries for building applications using mrgingham

%prep
%setup -q

%build
make %{?_smp_mflags} all

%install
rm -rf $RPM_BUILD_ROOT
%make_install

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%files
%doc
%{_libdir}/*.so.*
%{_bindir}/*
%{_mandir}/*
%{python3_sitelib}/*

%files devel
%{_libdir}/*.so
%{_includedir}/*
