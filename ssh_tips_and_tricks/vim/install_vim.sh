#!/bin/bash

# Clone the vim repository to install by hand, than you go inside the directory and run this script

# cop
setup python v3_9_2 # not sure if necessary
./configure --with-features=huge --prefix=$HOME/.local \
    --enable-multibyte \
        --enable-rubyinterp=yes \
        --enable-python3interp=yes \
        --with-python3-command=$PYTHON_VER \
        --with-python3-config-dir=$(python3-config --configdir) \
        --enable-perlinterp=yes \
        --enable-gui=gtk2 \
        --enable-cscope \

make

make install


# ADD $HOME/.local to your path on your .bashrc
