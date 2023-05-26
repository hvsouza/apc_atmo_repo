# How to use 

For the plugins set on the config file to work, clone the plugin menager of tmux as following:

``` sh
git clone https://github.com/tmux-plugins/tpm ~/.tmux/plugins/tpm
```


Open a tmux session and it should be all fine.

If you want _italic_ or **bold**, you need to setup a screen as following:

1. install the new screen `tic xterm-256color-italic.terminfo`
2. set it (put on your .bashrc): `export TERM=xterm-256color-italic`

This should work hehe

