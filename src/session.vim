let SessionLoad = 1
if &cp | set nocp | endif
nnoremap <silent>  :TlistToggle
nnoremap <silent>  :nohl
nnoremap <silent>  :set formatoptions+=at
nnoremap <silent> \b :CommandTBuffer
nnoremap <silent> \t :CommandT
let s:cpo_save=&cpo
set cpo&vim
nmap gx <Plug>NetrwBrowseX
nnoremap td :tabclose
nnoremap tm :tabm 
nnoremap tn :tabnext 
nnoremap tt :tabedit 
nnoremap tl :tablast
nnoremap tk :tabprev
nnoremap tj :tabnext
nnoremap th :tabfirst
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#NetrwBrowseX(expand("<cWORD>"),0)
let &cpo=s:cpo_save
unlet s:cpo_save
set autoindent
set background=dark
set backspace=indent,eol,start
set cscopeprg=/usr/bin/cscope
set cscopetag
set cscopeverbose
set expandtab
set fileencodings=ucs-bom,utf-8,latin1
set guicursor=n-v-c:block,o:hor50,i-ci:hor15,r-cr:hor30,sm:block,a:blinkon0
set helplang=en
set hlsearch
set ruler
set runtimepath=~/.vim,~/.vim/bundle/Vundle.vim,~/.vim/bundle/molokai,~/.vim/bundle/tabular,~/.vim/bundle/vim-fugitive,~/.vim/bundle/L9,~/.vim/bundle/command-t,/usr/share/vim/vimfiles,/usr/share/vim/vim72,/usr/share/vim/vimfiles/after,~/.vim/after,~/.vim/bundle/Vundle.vim,~/.vim/bundle/Vundle.vim/after,~/.vim/bundle/molokai/after,~/.vim/bundle/tabular/after,~/.vim/bundle/vim-fugitive/after,~/.vim/bundle/L9/after,~/.vim/bundle/command-t/after
set shiftwidth=2
set tabstop=2
set window=60
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Repos/openmc-mic-mpi/src
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +0 eigenvalue.F90
badd +137 initialize.F90
badd +0 tally_initialize.F90
badd +0 ~/Repos/openmc-base-perf/src/initialize.F90
args eigenvalue.F90
edit eigenvalue.F90
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
setlocal keymap=
setlocal noarabic
setlocal autoindent
setlocal nobinary
setlocal bufhidden=
setlocal buflisted
setlocal buftype=
setlocal nocindent
setlocal cinkeys=0{,0},0),:,0#,!^F,o,O,e
setlocal cinoptions=
setlocal cinwords=if,else,while,do,for,switch
setlocal comments=:!
setlocal commentstring=!%s
setlocal complete=.,w,b,u,t,i
setlocal completefunc=
setlocal nocopyindent
setlocal nocursorcolumn
setlocal nocursorline
setlocal define=
setlocal dictionary=
setlocal nodiff
setlocal equalprg=
setlocal errorformat=
setlocal expandtab
if &filetype != 'fortran'
setlocal filetype=fortran
endif
setlocal foldcolumn=0
setlocal foldenable
setlocal foldexpr=0
setlocal foldignore=#
setlocal foldlevel=2
setlocal foldmarker={{{,}}}
set foldmethod=syntax
setlocal foldmethod=syntax
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldtext=foldtext()
setlocal formatexpr=
setlocal formatoptions=tcql
setlocal formatlistpat=^\\s*\\d\\+[\\]:.)}\\t\ ]\\s*
setlocal grepprg=
setlocal iminsert=0
setlocal imsearch=0
setlocal include=^\\c#\\=\\s*include\\s\\+
setlocal includeexpr=
setlocal indentexpr=FortranGetFreeIndent()
setlocal indentkeys=0{,0},:,0#,!^F,o,O,e,=~end,=~case,=~if,=~else,=~do,=~where,=~elsewhere,=~select,=~endif,=~enddo,=~endwhere,=~endselect,=~type,=~interface
setlocal noinfercase
setlocal iskeyword=@,48-57,_,192-255
setlocal keywordprg=
setlocal nolinebreak
setlocal nolisp
setlocal nolist
setlocal makeprg=
setlocal matchpairs=(:),{:},[:]
setlocal modeline
setlocal modifiable
setlocal nrformats=octal,hex
setlocal nonumber
setlocal numberwidth=4
setlocal omnifunc=
setlocal path=
setlocal nopreserveindent
setlocal nopreviewwindow
setlocal quoteescape=\\
setlocal noreadonly
setlocal norightleft
setlocal rightleftcmd=search
setlocal noscrollbind
setlocal shiftwidth=2
setlocal noshortname
setlocal nosmartindent
setlocal softtabstop=0
setlocal nospell
setlocal spellcapcheck=[.?!]\\_[\\])'\"\	\ ]\\+
setlocal spellfile=
setlocal spelllang=en
setlocal statusline=
setlocal suffixesadd=.f95,.f90,.for,.f,.F,.f77,.ftn,.fpp
setlocal swapfile
setlocal synmaxcol=3000
if &syntax != 'fortran'
setlocal syntax=fortran
endif
setlocal tabstop=2
setlocal tags=~/Repos/openmc-mic-mpi/.git/fortran.tags,~/Repos/openmc-mic-mpi/.git/tags,./tags,./TAGS,tags,TAGS
setlocal textwidth=80
setlocal thesaurus=
setlocal nowinfixheight
setlocal nowinfixwidth
setlocal wrap
setlocal wrapmargin=0
1
normal zo
39
normal zo
271
normal zo
861
normal zo
1
normal zo
let s:l = 75 - ((28 * winheight(0) + 29) / 59)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
75
normal! 032l
tabedit initialize.F90
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 118 + 118) / 236)
exe 'vert 2resize ' . ((&columns * 117 + 118) / 236)
argglobal
setlocal keymap=
setlocal noarabic
setlocal autoindent
setlocal nobinary
setlocal bufhidden=
setlocal buflisted
setlocal buftype=
setlocal nocindent
setlocal cinkeys=0{,0},0),:,0#,!^F,o,O,e
setlocal cinoptions=
setlocal cinwords=if,else,while,do,for,switch
setlocal comments=:!
setlocal commentstring=!%s
setlocal complete=.,w,b,u,t,i
setlocal completefunc=
setlocal nocopyindent
setlocal nocursorcolumn
setlocal nocursorline
setlocal define=
setlocal dictionary=
setlocal nodiff
setlocal equalprg=
setlocal errorformat=
setlocal expandtab
if &filetype != 'fortran'
setlocal filetype=fortran
endif
setlocal foldcolumn=0
setlocal foldenable
setlocal foldexpr=0
setlocal foldignore=#
setlocal foldlevel=2
setlocal foldmarker={{{,}}}
set foldmethod=syntax
setlocal foldmethod=syntax
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldtext=foldtext()
setlocal formatexpr=
setlocal formatoptions=tcql
setlocal formatlistpat=^\\s*\\d\\+[\\]:.)}\\t\ ]\\s*
setlocal grepprg=
setlocal iminsert=0
setlocal imsearch=0
setlocal include=^\\c#\\=\\s*include\\s\\+
setlocal includeexpr=
setlocal indentexpr=FortranGetFreeIndent()
setlocal indentkeys=0{,0},:,0#,!^F,o,O,e,=~end,=~case,=~if,=~else,=~do,=~where,=~elsewhere,=~select,=~endif,=~enddo,=~endwhere,=~endselect,=~type,=~interface
setlocal noinfercase
setlocal iskeyword=@,48-57,_,192-255
setlocal keywordprg=
setlocal nolinebreak
setlocal nolisp
setlocal nolist
setlocal makeprg=
setlocal matchpairs=(:),{:},[:]
setlocal modeline
setlocal modifiable
setlocal nrformats=octal,hex
setlocal nonumber
setlocal numberwidth=4
setlocal omnifunc=
setlocal path=
setlocal nopreserveindent
setlocal nopreviewwindow
setlocal quoteescape=\\
setlocal noreadonly
setlocal norightleft
setlocal rightleftcmd=search
setlocal noscrollbind
setlocal shiftwidth=2
setlocal noshortname
setlocal nosmartindent
setlocal softtabstop=0
setlocal nospell
setlocal spellcapcheck=[.?!]\\_[\\])'\"\	\ ]\\+
setlocal spellfile=
setlocal spelllang=en
setlocal statusline=
setlocal suffixesadd=.f95,.f90,.for,.f,.F,.f77,.ftn,.fpp
setlocal swapfile
setlocal synmaxcol=3000
if &syntax != 'fortran'
setlocal syntax=fortran
endif
setlocal tabstop=2
setlocal tags=~/Repos/openmc-mic-mpi/.git/fortran.tags,~/Repos/openmc-mic-mpi/.git/tags,./tags,./TAGS,tags,TAGS
setlocal textwidth=80
setlocal thesaurus=
setlocal nowinfixheight
setlocal nowinfixwidth
setlocal wrap
setlocal wrapmargin=0
1
normal zo
49
normal zo
193
normal zo
325
normal zo
502
normal zo
578
normal zo
842
normal zo
882
normal zo
967
normal zo
1
normal zo
let s:l = 138 - ((29 * winheight(0) + 29) / 58)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
138
normal! 0
wincmd w
argglobal
edit ~/Repos/openmc-base-perf/src/initialize.F90
setlocal keymap=
setlocal noarabic
setlocal autoindent
setlocal nobinary
setlocal bufhidden=
setlocal buflisted
setlocal buftype=
setlocal nocindent
setlocal cinkeys=0{,0},0),:,0#,!^F,o,O,e
setlocal cinoptions=
setlocal cinwords=if,else,while,do,for,switch
setlocal comments=:!
setlocal commentstring=!%s
setlocal complete=.,w,b,u,t,i
setlocal completefunc=
setlocal nocopyindent
setlocal nocursorcolumn
setlocal nocursorline
setlocal define=
setlocal dictionary=
setlocal nodiff
setlocal equalprg=
setlocal errorformat=
setlocal expandtab
if &filetype != 'fortran'
setlocal filetype=fortran
endif
setlocal foldcolumn=0
setlocal foldenable
setlocal foldexpr=0
setlocal foldignore=#
setlocal foldlevel=2
setlocal foldmarker={{{,}}}
set foldmethod=syntax
setlocal foldmethod=syntax
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldtext=foldtext()
setlocal formatexpr=
setlocal formatoptions=tcql
setlocal formatlistpat=^\\s*\\d\\+[\\]:.)}\\t\ ]\\s*
setlocal grepprg=
setlocal iminsert=0
setlocal imsearch=0
setlocal include=^\\c#\\=\\s*include\\s\\+
setlocal includeexpr=
setlocal indentexpr=FortranGetFreeIndent()
setlocal indentkeys=0{,0},:,0#,!^F,o,O,e,=~end,=~case,=~if,=~else,=~do,=~where,=~elsewhere,=~select,=~endif,=~enddo,=~endwhere,=~endselect,=~type,=~interface
setlocal noinfercase
setlocal iskeyword=@,48-57,_,192-255
setlocal keywordprg=
setlocal nolinebreak
setlocal nolisp
setlocal nolist
setlocal makeprg=
setlocal matchpairs=(:),{:},[:]
setlocal modeline
setlocal modifiable
setlocal nrformats=octal,hex
setlocal nonumber
setlocal numberwidth=4
setlocal omnifunc=
setlocal path=
setlocal nopreserveindent
setlocal nopreviewwindow
setlocal quoteescape=\\
setlocal noreadonly
setlocal norightleft
setlocal rightleftcmd=search
setlocal noscrollbind
setlocal shiftwidth=2
setlocal noshortname
setlocal nosmartindent
setlocal softtabstop=0
setlocal nospell
setlocal spellcapcheck=[.?!]\\_[\\])'\"\	\ ]\\+
setlocal spellfile=
setlocal spelllang=en
setlocal statusline=
setlocal suffixesadd=.f95,.f90,.for,.f,.F,.f77,.ftn,.fpp
setlocal swapfile
setlocal synmaxcol=3000
if &syntax != 'fortran'
setlocal syntax=fortran
endif
setlocal tabstop=2
setlocal tags=~/Repos/openmc-base-perf/.git/fortran.tags,~/Repos/openmc-base-perf/.git/tags,./tags,./TAGS,tags,TAGS
setlocal textwidth=80
setlocal thesaurus=
setlocal nowinfixheight
setlocal nowinfixwidth
setlocal wrap
setlocal wrapmargin=0
1
normal zo
49
normal zo
822
normal zo
1
normal zo
let s:l = 818 - ((0 * winheight(0) + 29) / 58)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
818
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 118 + 118) / 236)
exe 'vert 2resize ' . ((&columns * 117 + 118) / 236)
tabedit tally_initialize.F90
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
setlocal keymap=
setlocal noarabic
setlocal autoindent
setlocal nobinary
setlocal bufhidden=
setlocal buflisted
setlocal buftype=
setlocal nocindent
setlocal cinkeys=0{,0},0),:,0#,!^F,o,O,e
setlocal cinoptions=
setlocal cinwords=if,else,while,do,for,switch
setlocal comments=:!
setlocal commentstring=!%s
setlocal complete=.,w,b,u,t,i
setlocal completefunc=
setlocal nocopyindent
setlocal nocursorcolumn
setlocal nocursorline
setlocal define=
setlocal dictionary=
setlocal nodiff
setlocal equalprg=
setlocal errorformat=
setlocal expandtab
if &filetype != 'fortran'
setlocal filetype=fortran
endif
setlocal foldcolumn=0
setlocal foldenable
setlocal foldexpr=0
setlocal foldignore=#
setlocal foldlevel=2
setlocal foldmarker={{{,}}}
set foldmethod=syntax
setlocal foldmethod=syntax
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldtext=foldtext()
setlocal formatexpr=
setlocal formatoptions=tcql
setlocal formatlistpat=^\\s*\\d\\+[\\]:.)}\\t\ ]\\s*
setlocal grepprg=
setlocal iminsert=0
setlocal imsearch=0
setlocal include=^\\c#\\=\\s*include\\s\\+
setlocal includeexpr=
setlocal indentexpr=FortranGetFreeIndent()
setlocal indentkeys=0{,0},:,0#,!^F,o,O,e,=~end,=~case,=~if,=~else,=~do,=~where,=~elsewhere,=~select,=~endif,=~enddo,=~endwhere,=~endselect,=~type,=~interface
setlocal noinfercase
setlocal iskeyword=@,48-57,_,192-255
setlocal keywordprg=
setlocal nolinebreak
setlocal nolisp
setlocal nolist
setlocal makeprg=
setlocal matchpairs=(:),{:},[:]
setlocal modeline
setlocal modifiable
setlocal nrformats=octal,hex
setlocal nonumber
setlocal numberwidth=4
setlocal omnifunc=
setlocal path=
setlocal nopreserveindent
setlocal nopreviewwindow
setlocal quoteescape=\\
setlocal noreadonly
setlocal norightleft
setlocal rightleftcmd=search
setlocal noscrollbind
setlocal shiftwidth=2
setlocal noshortname
setlocal nosmartindent
setlocal softtabstop=0
setlocal nospell
setlocal spellcapcheck=[.?!]\\_[\\])'\"\	\ ]\\+
setlocal spellfile=
setlocal spelllang=en
setlocal statusline=
setlocal suffixesadd=.f95,.f90,.for,.f,.F,.f77,.ftn,.fpp
setlocal swapfile
setlocal synmaxcol=3000
if &syntax != 'fortran'
setlocal syntax=fortran
endif
setlocal tabstop=2
setlocal tags=~/Repos/openmc-mic-mpi/.git/fortran.tags,~/Repos/openmc-mic-mpi/.git/tags,./tags,./TAGS,tags,TAGS
setlocal textwidth=80
setlocal thesaurus=
setlocal nowinfixheight
setlocal nowinfixwidth
setlocal wrap
setlocal wrapmargin=0
1
normal zo
20
normal zo
1
normal zo
let s:l = 88 - ((1 * winheight(0) + 21) / 42)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
88
normal! 013l
tabnext 2
if exists('s:wipebuf')
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToO
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
