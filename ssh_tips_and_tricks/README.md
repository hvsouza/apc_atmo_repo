# Some helpful tips or scripts for sshing

  * In this folder you will find to install vim and tmux, as much as examples of configuration file for your `.bashrc`.  
  (NOTE: you don't need to install vim or tmux, but if you want to customize, installing newer versions is quite mandatory)

## FERMILAB

* SSH on GPVM machines:
  * The `January 2023 LArSoft/Art tutorial` will give you all information you need.
  * For setting up a VNC connection, check out [this](https://wiki.dunescience.org/wiki/DUNE_Computing/Using_VNC_Connections_on_the_dunegpvms) page.

## CERN
* SSH on lxplus:
  * ssh -XY your_user@lxplus.cern.ch
  * Setting up VNC can be check [here](https://homepages.uc.edu/~schreihf/uchenry/post/vnc-to-lxplus/).



## General configuration

* You will find in repo an example of configuration, you can try copying it into your `.bash_profile` file.  
Please, keep in mind the following (from [this](https://linuxize.com/post/bashrc-vs-bash-profile/) site):

    > `.bash_profile` is read and executed when Bash is invoked as an interactive login shell, while `.bashrc` is executed for an interactive non-login shell.  
    > Use `.bash_profile` to run commands that should run only once, such as customizing the `$PATH` environment variable. 

    > Put the commands that should run every time you launch a new shell in the `.bashrc` file. This include your aliases and functions , custom prompts, history customization, and so on.

    > Typically, `~/.bash_profile` contains lines like below that source the `.bashrc` file. This means each time you log in to the terminal, both files are read and executed.

    > ``` sh
    > if [ -f ~/.bashrc ]; then
    >   . ~/.bashrc
    > fi
    > ```
  
## Vim configuration

* You cannot simply copy and paste my vim configuration, follow this tutorial first: <https://www.linode.com/docs/guides/introduction-to-vim-customization/>.
* After this, use my `.vimrc` and `.vimrc.plug` as reference to create your own. 

## Tmux
* The script might be outdated. You can use this as reference to change it: <https://frankenthal.dev/post/tmux-on-lxplus_pt0/>
* I believe this is easier, keep in mind that my configurations are optmize for vim keybindings.
* I use a different display for tmux that allows for italic and bold characteres. Make sure you disable it if you use my configuration

