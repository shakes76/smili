# Maintainer: shakes

pkgname=smili
pkgver=2.0
pkgrel=1
pkgdesc="SMILI and sMILX provides easy medical image processing and scientific visualisation."
arch=(x86_64)
url="https://github.com/shakes76/smili"
license=(BSD-3-Clause)
depends=(
  vtk                     # visualization
  qt6-base                # libvtkGUISupportQt.so etc. (5 direct libs, 6 total libs)
  qt6-declarative         # libvtkGUISupportQtQuick.so (1 direct lib, 1 total lib)
  itk                     # image processing
  zlib                    # 5 direct libs, 206 total libs
  # common data libraries
  expat                   # 1 direct lib, 103 total libs
  jsoncpp                 # 7 direct libs, 26 total libs
  libxml2                 # 3 direct libs, 13 total libs
)
makedepends=(
  # build system
  cmake
  ninja
  # graphical toolkits
  qt6-tools
  eigen
)
optdepends=(
  # additional tools not listed in makedepends
  'graphviz: drawing tools'
  'libglvnd: OpenGL rendering'  # checked at runtime rather than compile-time
  # bindings
  'java-runtime=11: java bindings'
  # graphical toolkits
  'qt6-declarative: QML plugin'
  # direct dependencies of "some" VTK libs/modules
  'libx11: rendering'
  'libxcursor: rendering'
  'fontconfig: rendering fonts with fontconfig support'
  'freetype2: rendering fonts'
  'gl2ps: rendering to PostScript, PDF, and SVG'
  'openvr: rendering for virtual reality'
  'openxr: rendering for virtual and augmented reality'
  'openimagedenoise: rendering with raytracing support'
  'ospray: rendering with raytracing support'
  'openmpi: OpenMPI support'
  'viskores: accelerators support'
  'ffmpeg: IO module'
  'imath: IO module'
  'liblas: IO module'
  'libogg: IO module'
  'libtheora: IO module'
)
# options=(staticlibs)
source=(
  $url/archive/refs/tags/v${pkgver}Full.zip
)
# sha256sums=('cee64b98d270ff7302daf1ef13458dff5d5ac1ecb45d47723835f7f7d562c989')
# b2sums=('6e7dab56c4f48d066ca44637f8839e9abc973d2831381ef5aec860aa1daa49ac9beca223e439c12dda2c33056b45395043f2edc54f55f9d169719b6c96499f40')

prepare() {
  cd ${pkgname^^}-${pkgver}Full
  _fast_float_version=$(pacman -Q fast_float | sed -e 's/.* //; s/-.*//g')
}

build() {
  cmake -B build -S ${pkgname^^}-${pkgver}Full -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr \
    -DCMAKE_INSTALL_LICENSEDIR=share/licenses/${pkgname} \
    -DCMAKE_SKIP_RPATH=OFF \
    -DCMAKE_SKIP_INSTALL_RPATH=ON \
    -DBUILD_SHARED_LIBS=ON \
    -Wno-dev
  cmake --build build
}

package() {
  DESTDIR="${pkgdir}" cmake --install build

  # Move the vtk.jar to the arch-specific location…
  install -dv "${pkgdir}"/usr/share/java/vtk
  mv -v "${pkgdir}"/usr/lib/java/vtk.jar "${pkgdir}"/usr/share/java/vtk
  # …and the libs to the proper place
  mv "${pkgdir}"/usr/lib/java/vtk-Linux-${CARCH}/*.so "${pkgdir}"/usr/lib/
  rmdir "${pkgdir}"/usr/lib/java/{vtk-Linux-${CARCH}/,}

  # byte-compile python modules since the CMake build does not do it
  local site_packages=$(python -c "import site; print(site.getsitepackages()[0])")
  python -m compileall -o 0 -o 1 -o 2 --hardlink-dupes -s "${pkgdir}" "${pkgdir}"${site_packages}

  # Remove third party CMake patching for older versions than ours
  rm -rv "${pkgdir}"/usr/lib/cmake/vtk/patches/3.*
}
