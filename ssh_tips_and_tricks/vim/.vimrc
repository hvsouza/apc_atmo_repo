" Set compatibility to Vim only.
set nocompatible

" Helps force plug-ins to load correctly when it is turned back on below.
filetype off

" Highlight cursor position
set cursorline


" For plug-ins to load correctly.
filetype plugin indent on

" Turn off modelines
set modelines=0

" Automatically wrap text that extends beyond the screen length.
set wrap
" Vim's auto indentation feature does not work properly with text copied from outside of Vim. Press the <F2> key to toggle paste mode on/off.
nnoremap <F2> :set invpaste paste?<CR>
imap <F2> <C-O>:set invpaste paste?<CR>
set pastetoggle=<F2>

" " Some personal remapings
" jk to exit insert mode
inoremap jk <Esc>
" nnoremap k gk
" nnoremap j gj
nnoremap <expr> j v:count == 0 ? 'gj' : "\<Esc>".v:count.'j'
nnoremap <expr> k v:count == 0 ? 'gk' : "\<Esc>".v:count.'k'
inoremap <c-e> <End>
inoremap <c-a> <Home>

" Uncomment below to set the max textwidth. Use a value corresponding to the width of your screen.
" set textwidth=79

" Hitting `o` will not add comment, hitting <Enter> in insert mode will
" Removed t and c (auto wrap text and comment)
" set formatoptions=qrn1j
" Doing like this avoid overwritting by synyax plugins
" autocmd BufNewFile,BufRead * set formatoptions-=c
autocmd FileType * set formatoptions=qrn1j

set smarttab
set tabstop=4
set shiftwidth=4
set softtabstop=4
" Spaces instead of tab for readl tab in insert mode: <Ctrl+V> Tab
set expandtab
" For indendint with `>` or `<`
set noshiftround

" Remap leader to space
nnoremap <SPACE> <Nop>
let mapleader=" "

" For the `best way to navigate files in vim` I remap buffer navigation
nnoremap <Leader>h :bprevious<CR>
nnoremap <Leader>l :bnext<CR>
nnoremap <Leader>k :bfirst<CR>
nnoremap <Leader>j :blast<CR>
nnoremap <Leader>bb :ls<CR>:b<Space>
" Go to previous buffer like doomemacs (not sure why two enters)
nnoremap <Leader>bl :ls<CR>:b<Space>#<CR><CR>
" Kill buffer like doomemacs
nnoremap <Leader>bk :bd<CR>
" To explore like doomemacs
nnoremap <Leader>. :e<Space>

" Display 5 lines above/below the cursor when scrolling with a mouse.
set scrolloff=5
" Fixes common backspace problems
set backspace=indent,eol,start

" Speed up scrolling in Vim
set ttyfast

" Status bar
set laststatus=2

" Display options
set showmode
set showcmd

" Highlight matching pairs of brackets. Use the '%' character to jump between them.
set matchpairs+=<:>

" Show line numbers
set number relativenumber

" Set status line display
" set statusline=%F%m%r%h%w\ [FORMAT=%{&ff}]\ [TYPE=%Y]\ [POS=%l,%v][%p%%]\ [BUFFER=%n]\ %{strftime('%c')}

" Encoding
set encoding=utf-8

" Highlight matching search patterns
set hlsearch

" Enable incremental search
set incsearch
" Include matching uppercase words with lowercase search term
set ignorecase
" Include only uppercase words with uppercase search term
set smartcase

" Uncomment to prevent non-normal modes showing in powerline and below powerline.
set noshowmode

" General settings from DT
set wildmenu                    " Display all matches when tab complete.
set hidden                      " Needed to keep multiple buffers open
set nobackup                    " No auto backups
set noswapfile                  " No swap
set t_Co=256                    " Set if term supports 256 colors.


" Store info from no more than 100 files at a time, 9999 lines of text, 100kb of data. Useful for copying large amounts of data between files.
set viminfo='100,<9999,s100

" Set .C file as cpp (my usual root script extension)
autocmd BufRead,BufNewFile *.C set filetype=cpp


" Vim was slow in fermilab due to clipboard.
" Running with VNC was fast, but without it was slow. This solves the proble:
" https://stackoverflow.com/questions/14635295/vim-takes-a-very-long-time-to-start-up
set clipboard^=exclude:.*

" Call the .vimrc.plug file
if filereadable(expand("~/.vimrc.plug"))
    source ~/.vimrc.plug
endif
" Remap for NERDTree
noremap <leader>n :NERDTreeFocus<CR>
noremap <C-n> :NERDTree<CR>
noremap <C-t> :NERDTreeToggle<CR>
let NERDTreeShowLineNumbers=1
let NERDTreeShowHidden=1
let NERDTreeMinimalUI = 1
let g:NERDTreeWinSize=38

" Allo mouse using inside vim, F3 to toogle it
set mouse=a
map <F3> <ESC>:exec &mouse!=""? "set mouse=" : "set mouse=nv"<CR>

" get rid of the annoying Ex mode
nnoremap Q <Nop>

" Disable numbers for better copying and paste
silent! map <F10> :set invnumber invrelativenumber <CR>

" " In case you want to copy to xclip
" vmap "+y :!xclip -f -sel clip <CR>
" map "+p :r!xclip -o -sel clip <CR>

" Smart dir behaviour that keeps current folder in the file i'm editting

set autochdir


" The lightline.vim theme
let g:lightline = {
            \ 'colorscheme': 'dracula',
            \ }

let g:dracula_italic = 1
let g:dracula_bold = 1
let g:dracula_colorterm = 1
colorscheme dracula


" Turn on syntax highlighting.
syntax enable


" Use a line cursor within insert mode and a block cursor everywhere else.
"
" Using iTerm2? Go-to preferences / profile / colors and disable the smart bar
" cursor color. Then pick a cursor and highlight color that matches your theme.
" That will ensure your cursor is always visible within insert mode.
"
" Reference chart of values:
"   Ps = 0  -> blinking block.
"   Ps = 1  -> blinking block (default).
"   Ps = 2  -> steady block.
"   Ps = 3  -> blinking underline.
"   Ps = 4  -> steady underline.
"   Ps = 5  -> blinking bar (xterm).
"   Ps = 6  -> steady bar (xterm).
let &t_SI = "\e[6 q"
let &t_EI = "\e[0 q"

" For italic to work
let &t_ZH="\e[3m"
let &t_ZR="\e[23m"

" Use <C-L> to clear the highlighting of :set hlsearch.
if maparg('<C-G>', 'n') ==# ''
  nnoremap <silent> <C-G> :nohlsearch<C-R>=has('diff')?'<Bar>diffupdate':''<CR><CR><C-L>
endif


" filenames like *.xml, *.html, *.xhtml, ...
" These are the file extensions where this plugin is enabled.
let g:closetag_filenames = '*.xml,*.html,*.xhtml,*.phtml'

" filetypes like xml, html, xhtml, ...
" These are the file types where this plugin is enabled.
let g:closetag_filetypes = 'xml,html,xhtml,phtml'

" integer value [0|1]
" This will make the list of non-closing tags case-sensitive (e.g. `<Link>` will be closed while `<link>` wont
let g:closetag_emptyTags_caseSensitive = 1

" Shortcut for closing tags, default is `>`
let g:closetag_shortcut = '>'

" Add > at current position without closing the current tag, default is ''
let g:closetag_close_shortcut = '<leader>>'



" Start interactive EasyAlign in visual mode (e.g. vipga)
xmap ga <Plug>(EasyAlign)

" Start interactive EasyAlign for a motion/text object (e.g. gaip)
nmap ga <Plug>(EasyAlign)


" For undo tree
nnoremap <C-x>u :UndotreeToggle<CR>

set undodir=$HOME/.undodir
set undofile
