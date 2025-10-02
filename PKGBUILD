# Maintainer: shakes

pkgname=smili
pkgver=2.1
pkgrel=1
pkgtype=Alpha
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
  libtiff
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
  $url/archive/refs/tags/v${pkgver}${pkgtype}.zip
)
sha256sums=('4c6e94b8cf3652f51301d8635221c4740a4f17f27e90eb5eea7dc2d30f9a96d9')

prepare() {
  cd ${pkgname}-${pkgver}${pkgtype}
  _fast_float_version=$(pacman -Q fast_float | sed -e 's/.* //; s/-.*//g')
}

build() {
  cmake -B build -S ${pkgname}-${pkgver}${pkgtype} -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr \
    -DCMAKE_SKIP_RPATH=OFF \
    -DCMAKE_SKIP_INSTALL_RPATH=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_PLUGINS=ON \
    -DBUILD_DICOM_PLUGIN=ON \
    -Wno-dev
  cmake --build build
}

package() {
  DESTDIR="${pkgdir}" cmake --install build
}
