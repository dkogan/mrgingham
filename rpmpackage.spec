# These xxx markers are to be replaced by git_build_rpm
Name:           xxx
Version:        xxx

Release:        1%{?dist}
Summary:        Dot (or grid) finder for camera calibrations

License:        LGPL-2.1+
URL:            https://github.com/dkogan/mrgingham/
Source0:        https://github.com/dkogan/mrgingham/archive/%{version}.tar.gz#/%{name}-%{version}.tar.gz

BuildRequires:  opencv-devel
BuildRequires:  boost-devel
BuildRequires:  chrpath


# to build the manpages I need to run 'mrgingham-observe-pixel-uncertainty
# --help'
BuildRequires:  python

%description
Library to find a grid of points; used for calibration routines

%package devel
Requires:       %{name}%{_isa} = %{version}-%{release}
Summary:        Development files for mrgingham

%description devel
Headers and libraries for building applications using mrgingham

%package tools
Requires:       %{name}%{_isa} = %{version}-%{release}
Summary:        Executable tools for mrgingham

%description tools
Executable tools for mrgingham

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

%files devel
%{_libdir}/*.so
%{_includedir}/*

%files tools
%{_bindir}/*
