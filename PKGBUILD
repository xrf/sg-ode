# Maintainer: Fei Yuan <yuan@nscl.msu.edu>
pkgname=sg-ode-git
pkgver=latest
pkgrel=1
pkgdesc="ODE solver by Shampine and Gordon"
arch=(i686 x86_64)
url=https://github.com/xrf/sg-ode
license=(LGPL)
makedepends=(git)
provides=(sg-ode)
conflicts=(sg-ode)
source=($pkgname::git://github.com/xrf/sg-ode)
sha256sums=(SKIP)

pkgver() {
    cd "$srcdir/$pkgname"
    s=`git 2>/dev/null describe --long --tags`
    if [ $? -eq 0 ]
    then
        printf '%s' "$s" | sed 's/^v//;s/\([^-]*-\)g/r\1/;s/-/./g'
    else
        printf 'r%s.%s' "`git rev-list --count HEAD`" \
               "`git rev-parse --short HEAD`"
    fi
}

build() {
    cd "$srcdir/$pkgname"
    make
}

package() {
    cd "$srcdir/$pkgname"
    make DESTDIR="$pkgdir" PREFIX=/usr install
    install -Dm644 LICENSE "$pkgdir/usr/share/licenses/$pkgname/LICENSE"
}
