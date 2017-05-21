cmake .. -G "Unix Makefiles" -DOPENSSL_ROOT_DIR=/usr/lib/openssl -DOPENSSL_LIBRARIES=/usr/lib/openssl/engines-1.0.0 -DOPENSSL_INCLUDE_DIR=/usr/include/openssl -DLIBUV_INCLUDE_DIR=/usr/local/include -DLIBUV_LIBRARY=/usr/local/lib

# $ export CC=/c/MinGW/bin/gcc.exe
# $ export CXX=/c/MinGW/bin/g++.exe
# $ cmake -G "MinGW Makefiles" ..

# I used this to update gcc and g++ in the windows Bash
# sudo add-apt-repository ppa:ubuntu-toolchain-r/test
# sudo apt-get update
# sudo apt-get install gcc-5 g++-5

# sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5