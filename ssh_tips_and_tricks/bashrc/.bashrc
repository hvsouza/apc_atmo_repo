# This is my personal config

alias dune_setup='source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh'
#export DUNETPC_VERSION=v09_22_02

alias gocodes='cd /dune/app/users/hsouza/'
alias godata='cd /dune/data/users/hsouza/'
alias ls='ls --color=auto'
alias ll='ls --color=auto -lhtr'
alias root='root -l'

alias myvim='vim --servername SERVER --remote'
alias startvim='vim --servername SERVER'

alias cppwd='pwd | xclip -sel c'
#to show the where I am
# export PS1="[\A]\w$ "
# PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\][\A]\w\[\033[00m\]$ '
export PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '

HISTSIZE=60000
HISTFILESIZE=60000

# for tmux italic
# Install:
#   tic screen-256color.terminfo
# Usage:
#   export TERM=screen-256color
# export TERM=screen-256color
export TERM=xterm-256color-italic

# alias toogledisplay='source ~/toogle_vnc.sh' # I use this, you will probably not :)

# this avoid window from looking weird
shopt -s checkwinsize

# prevent Ctrl-s from freezing terminal
stty -ixon

# for sarching forward with ctrl+t
bind "\C-t":forward-search-history


# Add .local/bin to PATH
PATH="${HOME}/.local/bin:${PATH}"
