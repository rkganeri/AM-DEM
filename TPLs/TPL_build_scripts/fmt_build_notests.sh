cmake \
  -DCMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-11 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="fmt_install" \
  -DFMT_TEST=Off\
  ..

make install
    
