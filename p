
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(/etc/vim)
fchdir() to previous dir
sourcing "/etc/vim/vimrc"
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 23: sourcing "/usr/share/vim/vim91/debian.vim"
finished sourcing /usr/share/vim/vim91/debian.vim
continuing in /etc/vim/vimrc
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 34: sourcing "/usr/share/vim/vim91/syntax/syntax.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 20: sourcing "/usr/share/vim/vim91/syntax/synload.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 22: sourcing "/usr/share/vim/vim91/syntax/syncolor.vim"
chdir(/usr/share/vim/vim91/colors/lists)
fchdir() to previous dir
line 57: sourcing "/usr/share/vim/vim91/colors/lists/default.vim"
finished sourcing /usr/share/vim/vim91/colors/lists/default.vim
continuing in /usr/share/vim/vim91/syntax/syncolor.vim
finished sourcing /usr/share/vim/vim91/syntax/syncolor.vim
continuing in /usr/share/vim/vim91/syntax/synload.vim
finished sourcing /usr/share/vim/vim91/syntax/synload.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/filetype.vim"
not found in 'runtimepath': "ftdetect/*.vim"
finished sourcing /usr/share/vim/vim91/filetype.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

Executing BufRead Autocommands for "*.json"
autocommand setf json

Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat | runtime! scripts.vim | endif

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat    && (expand("<amatch>") =~# '\.conf$'^I|| getline(1) =~ '^#' || getline(2) =~ '^#'^I|| getline(3) =~ '^#' || getline(4) =~ '^#'^I|| getline(5) =~ '^#') |   setf FALLBACK conf | endif

finished sourcing /usr/share/vim/vim91/syntax/syntax.vim
continuing in /etc/vim/vimrc
finished sourcing /etc/vim/vimrc
chdir(/home/upc_gcp)
fchdir() to previous dir
sourcing "$HOME/.vimrc"
chdir(/home/upc_gcp/.vim/autoload)
fchdir() to previous dir
line 2: sourcing "/home/upc_gcp/.vim/autoload/plug.vim"
finished sourcing /home/upc_gcp/.vim/autoload/plug.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 14: sourcing "/usr/share/vim/vim91/ftoff.vim"
finished sourcing /usr/share/vim/vim91/ftoff.vim
continuing in plug#end
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 88: sourcing "/usr/share/vim/vim91/filetype.vim"
chdir(/home/upc_gcp/.vim/plugged/vimtex/ftdetect)
fchdir() to previous dir
line 2918: sourcing "/home/upc_gcp/.vim/plugged/vimtex/ftdetect/cls.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/ftdetect/cls.vim
continuing in /usr/share/vim/vim91/filetype.vim
chdir(/home/upc_gcp/.vim/plugged/vimtex/ftdetect)
fchdir() to previous dir
line 2918: sourcing "/home/upc_gcp/.vim/plugged/vimtex/ftdetect/tex.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/ftdetect/tex.vim
continuing in /usr/share/vim/vim91/filetype.vim
chdir(/home/upc_gcp/.vim/plugged/vimtex/ftdetect)
fchdir() to previous dir
line 2918: sourcing "/home/upc_gcp/.vim/plugged/vimtex/ftdetect/tikz.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/ftdetect/tikz.vim
continuing in /usr/share/vim/vim91/filetype.vim
finished sourcing /usr/share/vim/vim91/filetype.vim
continuing in plug#end
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 88: sourcing "/usr/share/vim/vim91/ftplugin.vim"
finished sourcing /usr/share/vim/vim91/ftplugin.vim
continuing in plug#end
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 88: sourcing "/usr/share/vim/vim91/indent.vim"
finished sourcing /usr/share/vim/vim91/indent.vim
continuing in plug#end
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 25: sourcing "/usr/share/vim/vim91/syntax/syntax.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 16: sourcing "/usr/share/vim/vim91/syntax/nosyntax.vim"
Executing BufEnter Autocommands for "*"
autocommand syn clear

autocommand if exists("b:current_syntax") | unlet b:current_syntax | endif

finished sourcing /usr/share/vim/vim91/syntax/nosyntax.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 20: sourcing "/usr/share/vim/vim91/syntax/synload.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 22: sourcing "/usr/share/vim/vim91/syntax/syncolor.vim"
finished sourcing /usr/share/vim/vim91/syntax/syncolor.vim
continuing in /usr/share/vim/vim91/syntax/synload.vim
finished sourcing /usr/share/vim/vim91/syntax/synload.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

finished sourcing /usr/share/vim/vim91/syntax/syntax.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/filetype.vim"
finished sourcing /usr/share/vim/vim91/filetype.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/ftplugin.vim"
finished sourcing /usr/share/vim/vim91/ftplugin.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/indent.vim"
finished sourcing /usr/share/vim/vim91/indent.vim
continuing in /home/upc_gcp/.vimrc
finished sourcing $HOME/.vimrc
chdir(/home/upc_gcp/.vim)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim/plugged/vimtex/plugin)
fchdir() to previous dir
sourcing "/home/upc_gcp/.vim/plugged/vimtex/plugin/vimtex.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/plugin/vimtex.vim
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/plugin)
fchdir() to previous dir
sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/plugin/coc.vim"
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 38: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/rpc.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/rpc.vim
continuing in /home/upc_gcp/.vim/plugged/coc.nvim/plugin/coc.vim
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 54: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/util.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/util.vim
continuing in coc#rpc#start_server
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 58: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/client.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/client.vim
continuing in coc#rpc#start_server
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 10: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/api.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/api.vim
continuing in <SNR>19_Enable
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/plugin/coc.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/getscriptPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/getscriptPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/gzip.vim"
finished sourcing /usr/share/vim/vim91/plugin/gzip.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/logiPat.vim"
finished sourcing /usr/share/vim/vim91/plugin/logiPat.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/manpager.vim"
finished sourcing /usr/share/vim/vim91/plugin/manpager.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/matchparen.vim"
finished sourcing /usr/share/vim/vim91/plugin/matchparen.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/netrwPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/netrwPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/rrhelper.vim"
finished sourcing /usr/share/vim/vim91/plugin/rrhelper.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/spellfile.vim"
finished sourcing /usr/share/vim/vim91/plugin/spellfile.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/tarPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/tarPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/tohtml.vim"
finished sourcing /usr/share/vim/vim91/plugin/tohtml.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/vimballPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/vimballPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/zipPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/zipPlugin.vim
chdir(/home/upc_gcp/.vim/pack/tpope/start)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim/pack/tpope/start/commentary/plugin)
fchdir() to previous dir
sourcing "/home/upc_gcp/.vim/pack/tpope/start/commentary/plugin/commentary.vim"
finished sourcing /home/upc_gcp/.vim/pack/tpope/start/commentary/plugin/commentary.vim
not found in 'runtimepath': "plugin/**/*.vim"
Reading viminfo file "/home/upc_gcp/.viminfo" info oldfiles
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
                        "./RunConfig/Case_1.json" 
"./RunConfig/Case_1.json" 233L, 3028B
Reading viminfo file "/home/upc_gcp/.viminfo" marks
Executing BufRead Autocommands for "*.json"
autocommand setf json

Executing FileType Autocommands for "*"
autocommand call LoadFTPlugin()

chdir(/usr/share/vim/vim91/ftplugin)
fchdir() to previous dir
line 18: sourcing "/usr/share/vim/vim91/ftplugin/json.vim"
finished sourcing /usr/share/vim/vim91/ftplugin/json.vim
continuing in <SNR>15_LoadFTPlugin
Executing FileType Autocommands for "*"
autocommand call s:LoadIndent()

chdir(/usr/share/vim/vim91/indent)
fchdir() to previous dir
line 14: sourcing "/usr/share/vim/vim91/indent/json.vim"
finished sourcing /usr/share/vim/vim91/indent/json.vim
continuing in <SNR>16_LoadIndent
Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

Executing FileType Autocommands for "*"
autocommand call s:Autocmd('FileType', expand('<amatch>'), +expand('<abuf>'))

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat | runtime! scripts.vim | endif

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat    && (expand("<amatch>") =~# '\.conf$'^I|| getline(1) =~ '^#' || getline(2) =~ '^#'^I|| getline(3) =~ '^#' || getline(4) =~ '^#'^I|| getline(5) =~ '^#') |   setf FALLBACK conf | endif

Executing BufRead Autocommands for "*"
autocommand call s:Autocmd('BufCreate', +expand('<abuf>'))

Executing BufWinEnter Autocommands for "*"
autocommand call s:Autocmd('BufWinEnter', +expand('<abuf>'), win_getid(), coc#window#visible_range(win_getid()))

chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 0: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/window.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/window.vim
continuing in BufWinEnter Autocommands for "*"
Executing BufWinEnter Autocommands for "*"
autocommand call s:Highlight_Matching_Pair()

Executing BufEnter Autocommands for "*"
autocommand call s:HandleBufEnter(+expand('<abuf>'))

Executing BufEnter Autocommands for "*"
autocommand sil call s:LocalBrowse(expand("<amatch>"))

Executing VimEnter Autocommands for "*"
autocommand call s:VimEnter()

chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 3: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/compat.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/compat.vim
continuing in <SNR>19_VimEnter
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 2: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/hlgroup.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/hlgroup.vim
continuing in <SNR>19_Highlight
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 8: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/color.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/color.vim
continuing in <SNR>41_to_hex_color
Executing VimEnter Autocommands for "*"
autocommand sil call s:VimEnter(expand("<amatch>"))

chdir(./RunConfig)
fchdir() to previous dir
Executing CursorMoved Autocommands for "*"
autocommand call s:Autocmd('CursorMoved', +expand('<abuf>'), [line('.'), col('.')])

Executing CursorMoved Autocommands for "*"
autocommand call s:Highlight_Matching_Pair()

chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/syntax/syncolor.vim"
finished sourcing /usr/share/vim/vim91/syntax/syncolor.vim
SpecialKey     xxx term=bold ctermfg=81 guifg=Cyan
EndOfBuffer    xxx links to NonText
NonText        xxx term=bold ctermfg=12 gui=bold guifg=Blue
Directory      xxx term=bold ctermfg=159 guifg=Cyan
ErrorMsg       xxx term=standout ctermfg=15 ctermbg=1 guifg=White guibg=Red
IncSearch      xxx term=reverse cterm=reverse gui=reverse
Search         xxx term=reverse ctermfg=0 ctermbg=11 guifg=Black guibg=Yellow
CurSearch      xxx links to Search
MoreMsg        xxx term=bold ctermfg=121 gui=bold guifg=SeaGreen
ModeMsg        xxx term=bold cterm=bold gui=bold
LineNr         xxx term=underline ctermfg=11 guifg=Yellow
LineNrAbove    xxx cleared
LineNrBelow    xxx cleared
CursorLineNr   xxx term=bold cterm=underline ctermfg=11 gui=bold guifg=Yellow
CursorLineSign xxx links to SignColumn
CursorLineFold xxx links to FoldColumn
Question       xxx term=standout ctermfg=121 gui=bold guifg=Green
StatusLine     xxx term=bold,reverse cterm=bold,reverse gui=bold,reverse
StatusLineNC   xxx term=reverse cterm=reverse gui=reverse
VertSplit      xxx term=reverse cterm=reverse gui=reverse
Title          xxx term=bold ctermfg=225 gui=bold guifg=Magenta
Visual         xxx term=reverse ctermbg=242 guibg=DarkGrey
VisualNOS      xxx cleared
WarningMsg     xxx term=standout ctermfg=224 guifg=Red
WildMenu       xxx term=standout ctermfg=0 ctermbg=11 guifg=Black guibg=Yellow
Folded         xxx term=standout ctermfg=14 ctermbg=242 guifg=Cyan guibg=DarkGrey
FoldColumn     xxx term=standout ctermfg=14 ctermbg=242 guifg=Cyan guibg=Grey
DiffAdd        xxx term=bold ctermbg=4 guibg=DarkBlue
DiffChange     xxx term=bold ctermbg=5 guibg=DarkMagenta
DiffDelete     xxx term=bold ctermfg=12 ctermbg=6 gui=bold guifg=Blue guibg=DarkCyan
DiffText       xxx term=reverse cterm=bold ctermbg=9 gui=bold guibg=Red
SignColumn     xxx term=standout ctermfg=14 ctermbg=242 guifg=Cyan guibg=Grey
Conceal        xxx ctermfg=7 ctermbg=242 guifg=LightGrey guibg=DarkGrey
SpellBad       xxx term=reverse ctermbg=9 gui=undercurl guisp=Red
SpellCap       xxx term=reverse ctermbg=12 gui=undercurl guisp=Blue
SpellRare      xxx term=reverse ctermbg=13 gui=undercurl guisp=Magenta
SpellLocal     xxx term=underline ctermbg=14 gui=undercurl guisp=Cyan
Pmenu          xxx ctermfg=0 ctermbg=13 guibg=Magenta
PmenuSel       xxx ctermfg=242 ctermbg=0 guibg=DarkGrey
PmenuKind      xxx links to Pmenu
PmenuKindSel   xxx links to PmenuSel
PmenuExtra     xxx links to Pmenu
PmenuExtraSel  xxx links to PmenuSel
PmenuSbar      xxx ctermbg=248 guibg=Grey
PmenuThumb     xxx ctermbg=15 guibg=White
TabLine        xxx term=underline cterm=underline ctermfg=15 ctermbg=242 gui=underline guibg=DarkGrey
TabLineSel     xxx term=bold cterm=bold gui=bold
TabLineFill    xxx term=reverse cterm=reverse gui=reverse
CursorColumn   xxx term=reverse ctermbg=242 guibg=Grey40
CursorLine     xxx term=underline cterm=underline guibg=Grey40
ColorColumn    xxx term=reverse ctermbg=1 guibg=DarkRed
QuickFixLine   xxx links to Search
StatusLineTerm xxx term=bold,reverse cterm=bold ctermfg=0 ctermbg=121 gui=bold guifg=bg guibg=LightGreen
StatusLineTermNC xxx term=reverse ctermfg=0 ctermbg=121 guifg=bg guibg=LightGreen
Normal         xxx cleared
MatchParen     xxx term=reverse ctermbg=6 guibg=DarkCyan
ToolbarLine    xxx term=underline ctermbg=242 guibg=Grey50
ToolbarButton  xxx cterm=bold ctermfg=0 ctermbg=7 gui=bold guifg=Black guibg=LightGrey
Comment        xxx term=bold ctermfg=14 guifg=#80a0ff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 35
Constant       xxx term=underline ctermfg=13 guifg=#ffa0a0
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 36
Special        xxx term=bold ctermfg=224 guifg=Orange
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 37
Identifier     xxx term=underline cterm=bold ctermfg=14 guifg=#40ffff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 38
Statement      xxx term=bold ctermfg=11 gui=bold guifg=#ffff60
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 39
PreProc        xxx term=underline ctermfg=81 guifg=#ff80ff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 40
Type           xxx term=underline ctermfg=121 gui=bold guifg=#60ff60
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 41
Underlined     xxx term=underline cterm=underline ctermfg=81 gui=underline guifg=#80a0ff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 42
Ignore         xxx ctermfg=0 guifg=bg
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 43
Added          xxx ctermfg=10 guifg=LimeGreen
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 44
Changed        xxx ctermfg=12 guifg=DodgerBlue
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 45
Removed        xxx ctermfg=9 guifg=Red
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 46
Error          xxx term=reverse ctermfg=15 ctermbg=9 guifg=White guibg=Red
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 63
Todo           xxx term=standout ctermfg=0 ctermbg=11 guifg=Blue guibg=Yellow
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 64
String         xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 68
Character      xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 69
Number         xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 70
Boolean        xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 71
Float          xxx links to Number
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 72
Function       xxx links to Identifier
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 73
Conditional    xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 74
Repeat         xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 75
Label          xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 76
Operator       xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 77
Keyword        xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 78
Exception      xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 79
Include        xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 80
Define         xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 81
Macro          xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 82
PreCondit      xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 83
StorageClass   xxx links to Type
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 84
Structure      xxx links to Type
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 85
Typedef        xxx links to Type
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 86
Tag            xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 87
SpecialChar    xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 88
Delimiter      xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 89
SpecialComment xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 90
Debug          xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 91
jsonNoise      xxx links to Noise
	Last set from /usr/share/vim/vim91/syntax/json.vim line 121
jsonKeyword    xxx links to Label
	Last set from /usr/share/vim/vim91/syntax/json.vim line 108
jsonKeywordMatch xxx cleared
jsonQuote      xxx links to Quote
	Last set from /usr/share/vim/vim91/syntax/json.vim line 120
jsonString     xxx links to String
	Last set from /usr/share/vim/vim91/syntax/json.vim line 101
jsonStringMatch xxx cleared
jsonEscape     xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/json.vim line 103
jsonStringSQError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 116
jsonNumber     xxx links to Number
	Last set from /usr/share/vim/vim91/syntax/json.vim line 104
jsonNoQuotesError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 117
jsonTripleQuotesError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 118
jsonNumError   xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 111
jsonCommentError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 112
jsonSemicolonError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 113
jsonTrailingCommaError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 114
jsonMissingCommaError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 115
jsonPadding    xxx links to Operator
	Last set from /usr/share/vim/vim91/syntax/json.vim line 100
jsonBoolean    xxx links to Boolean
	Last set from /usr/share/vim/vim91/syntax/json.vim line 107
jsonNull       xxx links to Function
	Last set from /usr/share/vim/vim91/syntax/json.vim line 106
jsonBraces     xxx links to Delimiter
	Last set from /usr/share/vim/vim91/syntax/json.vim line 105
jsonFold       xxx cleared
jsonTest       xxx links to Label
	Last set from /usr/share/vim/vim91/syntax/json.vim line 102
Quote          xxx cleared
Noise          xxx cleared
CocSelectedText xxx ctermfg=9 guifg=#fb4934
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 459
CocCodeLens    xxx ctermfg=248 guifg=#999999
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 460
CocUnderline   xxx term=underline cterm=underline gui=underline guisp=#ebdbb2
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 461
CocBold        xxx term=bold cterm=bold gui=bold
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 462
CocItalic      xxx term=italic cterm=italic gui=italic
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 463
CocStrikeThrough xxx term=strikethrough cterm=strikethrough gui=strikethrough
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 464
CocMarkdownLink xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 465
CocDisabled    xxx ctermfg=248 guifg=#999999
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 466
CocSearch      xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 467
CocLink        xxx term=underline cterm=underline gui=underline guisp=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 468
CocFloatActive xxx links to CocSearch
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 469
CocFadeOut     xxx links to Conceal
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 470
CocMarkdownCode xxx links to markdownCode
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 471
markdownCode   xxx cleared
CocMarkdownHeader xxx links to markdownH1
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 472
markdownH1     xxx cleared
CocDeprecatedHighlight xxx links to CocStrikeThrough
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 473
CocUnusedHighlight xxx links to CocFadeOut
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 474
CocListSearch  xxx links to CocSearch
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 475
CocListMode    xxx links to ModeMsg
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 476
CocListPath    xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 477
CocHighlightText xxx links to CursorColumn
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 478
CocHoverRange  xxx links to Search
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 479
CocCursorRange xxx links to Search
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 480
CocLinkedEditing xxx links to CocCursorRange
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 481
CocHighlightRead xxx links to CocHighlightText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 482
CocHighlightWrite xxx links to CocHighlightText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 483
CocNotificationProgress xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 485
CocNotificationButton xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 486
CocNotificationKey xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 487
CocNotificationError xxx links to CocErrorFloat
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 488
CocErrorFloat  xxx ctermfg=9 ctermbg=253 guifg=#ff0000 guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocNotificationWarning xxx links to CocWarningFloat
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 489
CocWarningFloat xxx ctermfg=130 ctermbg=253 guifg=#ff922b guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocNotificationInfo xxx links to CocInfoFloat
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 490
CocInfoFloat   xxx ctermfg=11 ctermbg=253 guifg=#fab005 guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocSnippetVisual xxx links to Visual
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 492
CocTreeTitle   xxx links to Title
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 494
CocTreeDescription xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 495
CocTreeOpenClose xxx links to CocBold
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 496
CocTreeSelected xxx links to CursorLine
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 497
CocSelectedRange xxx links to CocHighlightText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 498
CocSymbolDefault xxx links to MoreMsg
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 500
CocPumSearch   xxx links to CocSearch
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 502
CocPumDetail   xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 503
CocPumMenu     xxx links to CocFloating
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 504
CocFloating    xxx ctermbg=253 guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 525
CocPumShortcut xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 505
CocPumDeprecated xxx links to CocStrikeThrough
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 506
CocVirtualText xxx ctermfg=12 guifg=#504945
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 507
CocPumVirtualText xxx links to CocVirtualText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 508
CocInputBoxVirtualText xxx links to CocVirtualText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 509
CocFloatDividingLine xxx links to CocVirtualText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 510
CocInlineVirtualText xxx ctermfg=244 guifg=#808080
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 512
CocInlineAnnotation xxx links to MoreMsg
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 516
CocListBlackBlack xxx guifg=#282828 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackBlue xxx guifg=#282828 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackGreen xxx guifg=#282828 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackGrey xxx guifg=#282828 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackWhite xxx guifg=#282828 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackCyan xxx guifg=#282828 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackYellow xxx guifg=#282828 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackMagenta xxx guifg=#282828 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackRed xxx guifg=#282828 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgBlack xxx ctermfg=0 guifg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgBlack xxx ctermbg=0 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListBlueBlack xxx guifg=#458588 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueBlue xxx guifg=#458588 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueGreen xxx guifg=#458588 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueGrey xxx guifg=#458588 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueWhite xxx guifg=#458588 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueCyan xxx guifg=#458588 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueYellow xxx guifg=#458588 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueMagenta xxx guifg=#458588 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueRed xxx guifg=#458588 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgBlue  xxx ctermfg=12 guifg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgBlue  xxx ctermbg=12 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListGreenBlack xxx guifg=#98971a guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenBlue xxx guifg=#98971a guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenGreen xxx guifg=#98971a guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenGrey xxx guifg=#98971a guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenWhite xxx guifg=#98971a guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenCyan xxx guifg=#98971a guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenYellow xxx guifg=#98971a guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenMagenta xxx guifg=#98971a guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenRed xxx guifg=#98971a guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgGreen xxx ctermfg=10 guifg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgGreen xxx ctermbg=10 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListGreyBlack xxx guifg=#928374 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyBlue xxx guifg=#928374 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyGreen xxx guifg=#928374 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyGrey xxx guifg=#928374 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyWhite xxx guifg=#928374 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyCyan xxx guifg=#928374 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyYellow xxx guifg=#928374 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyMagenta xxx guifg=#928374 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyRed xxx guifg=#928374 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgGrey  xxx ctermfg=248 guifg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgGrey  xxx ctermbg=248 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListWhiteBlack xxx guifg=#a89984 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteBlue xxx guifg=#a89984 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteGreen xxx guifg=#a89984 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteGrey xxx guifg=#a89984 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteWhite xxx guifg=#a89984 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteCyan xxx guifg=#a89984 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteYellow xxx guifg=#a89984 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteMagenta xxx guifg=#a89984 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteRed xxx guifg=#a89984 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgWhite xxx ctermfg=15 guifg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgWhite xxx ctermbg=15 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListCyanBlack xxx guifg=#689d6a guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanBlue xxx guifg=#689d6a guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanGreen xxx guifg=#689d6a guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanGrey xxx guifg=#689d6a guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanWhite xxx guifg=#689d6a guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanCyan xxx guifg=#689d6a guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanYellow xxx guifg=#689d6a guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanMagenta xxx guifg=#689d6a guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanRed xxx guifg=#689d6a guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgCyan  xxx ctermfg=14 guifg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgCyan  xxx ctermbg=14 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListYellowBlack xxx guifg=#d79921 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowBlue xxx guifg=#d79921 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowGreen xxx guifg=#d79921 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowGrey xxx guifg=#d79921 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowWhite xxx guifg=#d79921 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowCyan xxx guifg=#d79921 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowYellow xxx guifg=#d79921 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowMagenta xxx guifg=#d79921 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowRed xxx guifg=#d79921 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgYellow xxx ctermfg=11 guifg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgYellow xxx ctermbg=11 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListMagentaBlack xxx guifg=#b16286 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaBlue xxx guifg=#b16286 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaGreen xxx guifg=#b16286 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaGrey xxx guifg=#b16286 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaWhite xxx guifg=#b16286 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaCyan xxx guifg=#b16286 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaYellow xxx guifg=#b16286 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaMagenta xxx guifg=#b16286 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaRed xxx guifg=#b16286 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgMagenta xxx ctermfg=13 guifg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgMagenta xxx ctermbg=13 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListRedBlack xxx guifg=#cc241d guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedBlue xxx guifg=#cc241d guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedGreen xxx guifg=#cc241d guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedGrey xxx guifg=#cc241d guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedWhite xxx guifg=#cc241d guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedCyan xxx guifg=#cc241d guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedYellow xxx guifg=#cc241d guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedMagenta xxx guifg=#cc241d guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedRed  xxx guifg=#cc241d guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgRed   xxx ctermfg=9 guifg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgRed   xxx ctermbg=9 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocMenuSel     xxx ctermbg=250 guibg=#c8c8c8
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 526
CocListLine    xxx ctermbg=254 guibg=#eaeaea
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 542
CocFloatThumb  xxx ctermbg=249 guibg=#b7b7b7
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 527
CocFloatSbar   xxx links to CocFloating
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 528
CocFloatBorder xxx links to CocFloating
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 539
CocErrorHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocErrorSign   xxx ctermfg=9 guifg=#ff0000
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocErrorVirtualText xxx ctermfg=9 guifg=#ff0000
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocWarningHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocWarningSign xxx ctermfg=130 guifg=#ff922b
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocWarningVirtualText xxx ctermfg=130 guifg=#ff922b
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocInfoHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocInfoSign    xxx ctermfg=11 guifg=#fab005
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocInfoVirtualText xxx ctermfg=11 guifg=#fab005
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocHintHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocHintSign    xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocHintVirtualText xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocHintFloat   xxx ctermfg=12 ctermbg=253 guifg=#15aabf guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocInlayHint   xxx ctermfg=12 ctermbg=248 guifg=#15aabf guibg=Grey
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocInlayHintParameter xxx links to CocInlayHint
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 583
CocInlayHintType xxx links to CocInlayHint
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 583
CocSemTypeMacro xxx links to Define
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeEnumMember xxx links to Constant
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeComment xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeOperator xxx links to Operator
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeProperty xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeClass xxx links to Special
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeStruct xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeRegexp xxx links to String
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeBoolean xxx links to Boolean
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeMethod xxx links to Function
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeNamespace xxx links to Include
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeModifier xxx links to StorageClass
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeNumber xxx links to Number
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeTypeParameter xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeKeyword xxx links to Keyword
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemModDeprecated xxx links to CocDeprecatedHighlight
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeFunction xxx links to Function
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeDecorator xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeEnum xxx links to Type
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeParameter xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeType xxx links to Type
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeString xxx links to String
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeVariable xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeInterface xxx links to Type
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeEvent xxx links to Keyword
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSymbolUnit  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolNumber xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolFunction xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolKey   xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolKeyword xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolReference xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolFolder xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolVariable xxx ctermfg=5 guifg=#6a5acd
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolNull  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolValue xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolConstant xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolText  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolModule xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolPackage xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolClass xxx ctermfg=5 guifg=#6a5acd
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolOperator xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolStruct xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolObject xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolMethod xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolArray xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolEnum  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolField xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolInterface xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolProperty xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolColor xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolFile  xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolEvent xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolTypeParameter xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolConstructor xxx ctermfg=5 guifg=#6a5acd
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolSnippet xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolBoolean xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolNamespace xxx ctermfg=5 guifg=#6a0dad
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolString xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolEnumMember xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSelectedLine xxx cleared
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 10: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/ui.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/ui.vim
continuing in coc#api#Call_function
Executing ModeChanged Autocommands for "*"
autocommand call s:Autocmd('ModeChanged', v:event)


E173: 23 more files to edit
Executing ModeChanged Autocommands for "*"
autocommand call s:Autocmd('ModeChanged', v:event)


chdir(./RunConfig)
fchdir() to previous dir
Executing BufWinLeave Autocommands for "*"
autocommand call s:Autocmd('BufWinLeave', +expand('<abuf>'), bufwinid(+expand('<abuf>')))

chdir(./RunConfig)
fchdir() to previous dir
Executing BufUnload Autocommands for "*"
autocommand call s:Autocmd('BufUnload', +expand('<abuf>'))

Executing VimLeavePre Autocommands for "*"
autocommand call s:VimLeavePre()

Writing viminfo file "/home/upc_gcp/.viminfo"
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(/etc/vim)
fchdir() to previous dir
sourcing "/etc/vim/vimrc"
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 23: sourcing "/usr/share/vim/vim91/debian.vim"
finished sourcing /usr/share/vim/vim91/debian.vim
continuing in /etc/vim/vimrc
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 34: sourcing "/usr/share/vim/vim91/syntax/syntax.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 20: sourcing "/usr/share/vim/vim91/syntax/synload.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 22: sourcing "/usr/share/vim/vim91/syntax/syncolor.vim"
chdir(/usr/share/vim/vim91/colors/lists)
fchdir() to previous dir
line 57: sourcing "/usr/share/vim/vim91/colors/lists/default.vim"
finished sourcing /usr/share/vim/vim91/colors/lists/default.vim
continuing in /usr/share/vim/vim91/syntax/syncolor.vim
finished sourcing /usr/share/vim/vim91/syntax/syncolor.vim
continuing in /usr/share/vim/vim91/syntax/synload.vim
finished sourcing /usr/share/vim/vim91/syntax/synload.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/filetype.vim"
not found in 'runtimepath': "ftdetect/*.vim"
finished sourcing /usr/share/vim/vim91/filetype.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

Executing BufRead Autocommands for "*.json"
autocommand setf json

Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat | runtime! scripts.vim | endif

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat    && (expand("<amatch>") =~# '\.conf$'^I|| getline(1) =~ '^#' || getline(2) =~ '^#'^I|| getline(3) =~ '^#' || getline(4) =~ '^#'^I|| getline(5) =~ '^#') |   setf FALLBACK conf | endif

finished sourcing /usr/share/vim/vim91/syntax/syntax.vim
continuing in /etc/vim/vimrc
finished sourcing /etc/vim/vimrc
chdir(/home/upc_gcp)
fchdir() to previous dir
sourcing "$HOME/.vimrc"
chdir(/home/upc_gcp/.vim/autoload)
fchdir() to previous dir
line 2: sourcing "/home/upc_gcp/.vim/autoload/plug.vim"
finished sourcing /home/upc_gcp/.vim/autoload/plug.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 14: sourcing "/usr/share/vim/vim91/ftoff.vim"
finished sourcing /usr/share/vim/vim91/ftoff.vim
continuing in plug#end
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 88: sourcing "/usr/share/vim/vim91/filetype.vim"
chdir(/home/upc_gcp/.vim/plugged/vimtex/ftdetect)
fchdir() to previous dir
line 2918: sourcing "/home/upc_gcp/.vim/plugged/vimtex/ftdetect/cls.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/ftdetect/cls.vim
continuing in /usr/share/vim/vim91/filetype.vim
chdir(/home/upc_gcp/.vim/plugged/vimtex/ftdetect)
fchdir() to previous dir
line 2918: sourcing "/home/upc_gcp/.vim/plugged/vimtex/ftdetect/tex.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/ftdetect/tex.vim
continuing in /usr/share/vim/vim91/filetype.vim
chdir(/home/upc_gcp/.vim/plugged/vimtex/ftdetect)
fchdir() to previous dir
line 2918: sourcing "/home/upc_gcp/.vim/plugged/vimtex/ftdetect/tikz.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/ftdetect/tikz.vim
continuing in /usr/share/vim/vim91/filetype.vim
finished sourcing /usr/share/vim/vim91/filetype.vim
continuing in plug#end
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 88: sourcing "/usr/share/vim/vim91/ftplugin.vim"
finished sourcing /usr/share/vim/vim91/ftplugin.vim
continuing in plug#end
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 88: sourcing "/usr/share/vim/vim91/indent.vim"
finished sourcing /usr/share/vim/vim91/indent.vim
continuing in plug#end
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 25: sourcing "/usr/share/vim/vim91/syntax/syntax.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 16: sourcing "/usr/share/vim/vim91/syntax/nosyntax.vim"
Executing BufEnter Autocommands for "*"
autocommand syn clear

autocommand if exists("b:current_syntax") | unlet b:current_syntax | endif

finished sourcing /usr/share/vim/vim91/syntax/nosyntax.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 20: sourcing "/usr/share/vim/vim91/syntax/synload.vim"
chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
line 22: sourcing "/usr/share/vim/vim91/syntax/syncolor.vim"
finished sourcing /usr/share/vim/vim91/syntax/syncolor.vim
continuing in /usr/share/vim/vim91/syntax/synload.vim
finished sourcing /usr/share/vim/vim91/syntax/synload.vim
continuing in /usr/share/vim/vim91/syntax/syntax.vim
Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

finished sourcing /usr/share/vim/vim91/syntax/syntax.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/filetype.vim"
finished sourcing /usr/share/vim/vim91/filetype.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/ftplugin.vim"
finished sourcing /usr/share/vim/vim91/ftplugin.vim
continuing in /home/upc_gcp/.vimrc
chdir(/usr/share/vim/vim91)
fchdir() to previous dir
line 26: sourcing "/usr/share/vim/vim91/indent.vim"
finished sourcing /usr/share/vim/vim91/indent.vim
continuing in /home/upc_gcp/.vimrc
finished sourcing $HOME/.vimrc
chdir(/home/upc_gcp/.vim)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim/plugged/vimtex/plugin)
fchdir() to previous dir
sourcing "/home/upc_gcp/.vim/plugged/vimtex/plugin/vimtex.vim"
finished sourcing /home/upc_gcp/.vim/plugged/vimtex/plugin/vimtex.vim
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/plugin)
fchdir() to previous dir
sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/plugin/coc.vim"
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 38: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/rpc.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/rpc.vim
continuing in /home/upc_gcp/.vim/plugged/coc.nvim/plugin/coc.vim
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 54: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/util.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/util.vim
continuing in coc#rpc#start_server
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 58: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/client.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/client.vim
continuing in coc#rpc#start_server
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 10: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/api.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/api.vim
continuing in <SNR>19_Enable
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/plugin/coc.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/getscriptPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/getscriptPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/gzip.vim"
finished sourcing /usr/share/vim/vim91/plugin/gzip.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/logiPat.vim"
finished sourcing /usr/share/vim/vim91/plugin/logiPat.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/manpager.vim"
finished sourcing /usr/share/vim/vim91/plugin/manpager.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/matchparen.vim"
finished sourcing /usr/share/vim/vim91/plugin/matchparen.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/netrwPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/netrwPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/rrhelper.vim"
finished sourcing /usr/share/vim/vim91/plugin/rrhelper.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/spellfile.vim"
finished sourcing /usr/share/vim/vim91/plugin/spellfile.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/tarPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/tarPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/tohtml.vim"
finished sourcing /usr/share/vim/vim91/plugin/tohtml.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/vimballPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/vimballPlugin.vim
chdir(/usr/share/vim/vim91/plugin)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/plugin/zipPlugin.vim"
finished sourcing /usr/share/vim/vim91/plugin/zipPlugin.vim
chdir(/home/upc_gcp/.vim/pack/tpope/start)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim/pack/tpope/start/commentary/plugin)
fchdir() to previous dir
sourcing "/home/upc_gcp/.vim/pack/tpope/start/commentary/plugin/commentary.vim"
finished sourcing /home/upc_gcp/.vim/pack/tpope/start/commentary/plugin/commentary.vim
not found in 'runtimepath': "plugin/**/*.vim"
Reading viminfo file "/home/upc_gcp/.viminfo" info oldfiles
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
                        "./RunConfig/Case_1.json" 
"./RunConfig/Case_1.json" 233L, 3028B
Reading viminfo file "/home/upc_gcp/.viminfo" marks
Executing BufRead Autocommands for "*.json"
autocommand setf json

Executing FileType Autocommands for "*"
autocommand call LoadFTPlugin()

chdir(/usr/share/vim/vim91/ftplugin)
fchdir() to previous dir
line 18: sourcing "/usr/share/vim/vim91/ftplugin/json.vim"
finished sourcing /usr/share/vim/vim91/ftplugin/json.vim
continuing in <SNR>15_LoadFTPlugin
Executing FileType Autocommands for "*"
autocommand call s:LoadIndent()

chdir(/usr/share/vim/vim91/indent)
fchdir() to previous dir
line 14: sourcing "/usr/share/vim/vim91/indent/json.vim"
finished sourcing /usr/share/vim/vim91/indent/json.vim
continuing in <SNR>16_LoadIndent
Executing FileType Autocommands for "*"
autocommand 0verbose exe "set syntax=" . expand("<amatch>")

Executing FileType Autocommands for "*"
autocommand call s:Autocmd('FileType', expand('<amatch>'), +expand('<abuf>'))

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat | runtime! scripts.vim | endif

Executing BufRead Autocommands for "*"
autocommand if !did_filetype() && expand("<amatch>") !~ g:ft_ignore_pat    && (expand("<amatch>") =~# '\.conf$'^I|| getline(1) =~ '^#' || getline(2) =~ '^#'^I|| getline(3) =~ '^#' || getline(4) =~ '^#'^I|| getline(5) =~ '^#') |   setf FALLBACK conf | endif

Executing BufRead Autocommands for "*"
autocommand call s:Autocmd('BufCreate', +expand('<abuf>'))

Executing BufWinEnter Autocommands for "*"
autocommand call s:Autocmd('BufWinEnter', +expand('<abuf>'), win_getid(), coc#window#visible_range(win_getid()))

chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 0: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/window.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/window.vim
continuing in BufWinEnter Autocommands for "*"
Executing BufWinEnter Autocommands for "*"
autocommand call s:Highlight_Matching_Pair()

Executing BufEnter Autocommands for "*"
autocommand call s:HandleBufEnter(+expand('<abuf>'))

Executing BufEnter Autocommands for "*"
autocommand sil call s:LocalBrowse(expand("<amatch>"))

Executing VimEnter Autocommands for "*"
autocommand call s:VimEnter()

chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 3: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/compat.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/compat.vim
continuing in <SNR>19_VimEnter
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 2: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/hlgroup.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/hlgroup.vim
continuing in <SNR>19_Highlight
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 8: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/color.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/color.vim
continuing in <SNR>41_to_hex_color
Executing VimEnter Autocommands for "*"
autocommand sil call s:VimEnter(expand("<amatch>"))

chdir(./RunConfig)
fchdir() to previous dir
Executing CursorMoved Autocommands for "*"
autocommand call s:Autocmd('CursorMoved', +expand('<abuf>'), [line('.'), col('.')])

Executing CursorMoved Autocommands for "*"
autocommand call s:Highlight_Matching_Pair()

chdir(/usr/share/vim/vim91/syntax)
fchdir() to previous dir
sourcing "/usr/share/vim/vim91/syntax/syncolor.vim"
finished sourcing /usr/share/vim/vim91/syntax/syncolor.vim
SpecialKey     xxx term=bold ctermfg=81 guifg=Cyan
EndOfBuffer    xxx links to NonText
NonText        xxx term=bold ctermfg=12 gui=bold guifg=Blue
Directory      xxx term=bold ctermfg=159 guifg=Cyan
ErrorMsg       xxx term=standout ctermfg=15 ctermbg=1 guifg=White guibg=Red
IncSearch      xxx term=reverse cterm=reverse gui=reverse
Search         xxx term=reverse ctermfg=0 ctermbg=11 guifg=Black guibg=Yellow
CurSearch      xxx links to Search
MoreMsg        xxx term=bold ctermfg=121 gui=bold guifg=SeaGreen
ModeMsg        xxx term=bold cterm=bold gui=bold
LineNr         xxx term=underline ctermfg=11 guifg=Yellow
LineNrAbove    xxx cleared
LineNrBelow    xxx cleared
CursorLineNr   xxx term=bold cterm=underline ctermfg=11 gui=bold guifg=Yellow
CursorLineSign xxx links to SignColumn
CursorLineFold xxx links to FoldColumn
Question       xxx term=standout ctermfg=121 gui=bold guifg=Green
StatusLine     xxx term=bold,reverse cterm=bold,reverse gui=bold,reverse
StatusLineNC   xxx term=reverse cterm=reverse gui=reverse
VertSplit      xxx term=reverse cterm=reverse gui=reverse
Title          xxx term=bold ctermfg=225 gui=bold guifg=Magenta
Visual         xxx term=reverse ctermbg=242 guibg=DarkGrey
VisualNOS      xxx cleared
WarningMsg     xxx term=standout ctermfg=224 guifg=Red
WildMenu       xxx term=standout ctermfg=0 ctermbg=11 guifg=Black guibg=Yellow
Folded         xxx term=standout ctermfg=14 ctermbg=242 guifg=Cyan guibg=DarkGrey
FoldColumn     xxx term=standout ctermfg=14 ctermbg=242 guifg=Cyan guibg=Grey
DiffAdd        xxx term=bold ctermbg=4 guibg=DarkBlue
DiffChange     xxx term=bold ctermbg=5 guibg=DarkMagenta
DiffDelete     xxx term=bold ctermfg=12 ctermbg=6 gui=bold guifg=Blue guibg=DarkCyan
DiffText       xxx term=reverse cterm=bold ctermbg=9 gui=bold guibg=Red
SignColumn     xxx term=standout ctermfg=14 ctermbg=242 guifg=Cyan guibg=Grey
Conceal        xxx ctermfg=7 ctermbg=242 guifg=LightGrey guibg=DarkGrey
SpellBad       xxx term=reverse ctermbg=9 gui=undercurl guisp=Red
SpellCap       xxx term=reverse ctermbg=12 gui=undercurl guisp=Blue
SpellRare      xxx term=reverse ctermbg=13 gui=undercurl guisp=Magenta
SpellLocal     xxx term=underline ctermbg=14 gui=undercurl guisp=Cyan
Pmenu          xxx ctermfg=0 ctermbg=13 guibg=Magenta
PmenuSel       xxx ctermfg=242 ctermbg=0 guibg=DarkGrey
PmenuKind      xxx links to Pmenu
PmenuKindSel   xxx links to PmenuSel
PmenuExtra     xxx links to Pmenu
PmenuExtraSel  xxx links to PmenuSel
PmenuSbar      xxx ctermbg=248 guibg=Grey
PmenuThumb     xxx ctermbg=15 guibg=White
TabLine        xxx term=underline cterm=underline ctermfg=15 ctermbg=242 gui=underline guibg=DarkGrey
TabLineSel     xxx term=bold cterm=bold gui=bold
TabLineFill    xxx term=reverse cterm=reverse gui=reverse
CursorColumn   xxx term=reverse ctermbg=242 guibg=Grey40
CursorLine     xxx term=underline cterm=underline guibg=Grey40
ColorColumn    xxx term=reverse ctermbg=1 guibg=DarkRed
QuickFixLine   xxx links to Search
StatusLineTerm xxx term=bold,reverse cterm=bold ctermfg=0 ctermbg=121 gui=bold guifg=bg guibg=LightGreen
StatusLineTermNC xxx term=reverse ctermfg=0 ctermbg=121 guifg=bg guibg=LightGreen
Normal         xxx cleared
MatchParen     xxx term=reverse ctermbg=6 guibg=DarkCyan
ToolbarLine    xxx term=underline ctermbg=242 guibg=Grey50
ToolbarButton  xxx cterm=bold ctermfg=0 ctermbg=7 gui=bold guifg=Black guibg=LightGrey
Comment        xxx term=bold ctermfg=14 guifg=#80a0ff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 35
Constant       xxx term=underline ctermfg=13 guifg=#ffa0a0
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 36
Special        xxx term=bold ctermfg=224 guifg=Orange
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 37
Identifier     xxx term=underline cterm=bold ctermfg=14 guifg=#40ffff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 38
Statement      xxx term=bold ctermfg=11 gui=bold guifg=#ffff60
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 39
PreProc        xxx term=underline ctermfg=81 guifg=#ff80ff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 40
Type           xxx term=underline ctermfg=121 gui=bold guifg=#60ff60
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 41
Underlined     xxx term=underline cterm=underline ctermfg=81 gui=underline guifg=#80a0ff
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 42
Ignore         xxx ctermfg=0 guifg=bg
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 43
Added          xxx ctermfg=10 guifg=LimeGreen
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 44
Changed        xxx ctermfg=12 guifg=DodgerBlue
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 45
Removed        xxx ctermfg=9 guifg=Red
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 46
Error          xxx term=reverse ctermfg=15 ctermbg=9 guifg=White guibg=Red
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 63
Todo           xxx term=standout ctermfg=0 ctermbg=11 guifg=Blue guibg=Yellow
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 64
String         xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 68
Character      xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 69
Number         xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 70
Boolean        xxx links to Constant
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 71
Float          xxx links to Number
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 72
Function       xxx links to Identifier
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 73
Conditional    xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 74
Repeat         xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 75
Label          xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 76
Operator       xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 77
Keyword        xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 78
Exception      xxx links to Statement
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 79
Include        xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 80
Define         xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 81
Macro          xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 82
PreCondit      xxx links to PreProc
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 83
StorageClass   xxx links to Type
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 84
Structure      xxx links to Type
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 85
Typedef        xxx links to Type
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 86
Tag            xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 87
SpecialChar    xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 88
Delimiter      xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 89
SpecialComment xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 90
Debug          xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/syncolor.vim line 91
jsonNoise      xxx links to Noise
	Last set from /usr/share/vim/vim91/syntax/json.vim line 121
jsonKeyword    xxx links to Label
	Last set from /usr/share/vim/vim91/syntax/json.vim line 108
jsonKeywordMatch xxx cleared
jsonQuote      xxx links to Quote
	Last set from /usr/share/vim/vim91/syntax/json.vim line 120
jsonString     xxx links to String
	Last set from /usr/share/vim/vim91/syntax/json.vim line 101
jsonStringMatch xxx cleared
jsonEscape     xxx links to Special
	Last set from /usr/share/vim/vim91/syntax/json.vim line 103
jsonStringSQError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 116
jsonNumber     xxx links to Number
	Last set from /usr/share/vim/vim91/syntax/json.vim line 104
jsonNoQuotesError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 117
jsonTripleQuotesError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 118
jsonNumError   xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 111
jsonCommentError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 112
jsonSemicolonError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 113
jsonTrailingCommaError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 114
jsonMissingCommaError xxx links to Error
	Last set from /usr/share/vim/vim91/syntax/json.vim line 115
jsonPadding    xxx links to Operator
	Last set from /usr/share/vim/vim91/syntax/json.vim line 100
jsonBoolean    xxx links to Boolean
	Last set from /usr/share/vim/vim91/syntax/json.vim line 107
jsonNull       xxx links to Function
	Last set from /usr/share/vim/vim91/syntax/json.vim line 106
jsonBraces     xxx links to Delimiter
	Last set from /usr/share/vim/vim91/syntax/json.vim line 105
jsonFold       xxx cleared
jsonTest       xxx links to Label
	Last set from /usr/share/vim/vim91/syntax/json.vim line 102
Quote          xxx cleared
Noise          xxx cleared
CocSelectedText xxx ctermfg=9 guifg=#fb4934
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 459
CocCodeLens    xxx ctermfg=248 guifg=#999999
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 460
CocUnderline   xxx term=underline cterm=underline gui=underline guisp=#ebdbb2
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 461
CocBold        xxx term=bold cterm=bold gui=bold
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 462
CocItalic      xxx term=italic cterm=italic gui=italic
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 463
CocStrikeThrough xxx term=strikethrough cterm=strikethrough gui=strikethrough
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 464
CocMarkdownLink xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 465
CocDisabled    xxx ctermfg=248 guifg=#999999
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 466
CocSearch      xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 467
CocLink        xxx term=underline cterm=underline gui=underline guisp=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 468
CocFloatActive xxx links to CocSearch
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 469
CocFadeOut     xxx links to Conceal
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 470
CocMarkdownCode xxx links to markdownCode
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 471
markdownCode   xxx cleared
CocMarkdownHeader xxx links to markdownH1
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 472
markdownH1     xxx cleared
CocDeprecatedHighlight xxx links to CocStrikeThrough
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 473
CocUnusedHighlight xxx links to CocFadeOut
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 474
CocListSearch  xxx links to CocSearch
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 475
CocListMode    xxx links to ModeMsg
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 476
CocListPath    xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 477
CocHighlightText xxx links to CursorColumn
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 478
CocHoverRange  xxx links to Search
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 479
CocCursorRange xxx links to Search
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 480
CocLinkedEditing xxx links to CocCursorRange
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 481
CocHighlightRead xxx links to CocHighlightText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 482
CocHighlightWrite xxx links to CocHighlightText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 483
CocNotificationProgress xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 485
CocNotificationButton xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 486
CocNotificationKey xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 487
CocNotificationError xxx links to CocErrorFloat
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 488
CocErrorFloat  xxx ctermfg=9 ctermbg=253 guifg=#ff0000 guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocNotificationWarning xxx links to CocWarningFloat
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 489
CocWarningFloat xxx ctermfg=130 ctermbg=253 guifg=#ff922b guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocNotificationInfo xxx links to CocInfoFloat
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 490
CocInfoFloat   xxx ctermfg=11 ctermbg=253 guifg=#fab005 guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocSnippetVisual xxx links to Visual
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 492
CocTreeTitle   xxx links to Title
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 494
CocTreeDescription xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 495
CocTreeOpenClose xxx links to CocBold
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 496
CocTreeSelected xxx links to CursorLine
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 497
CocSelectedRange xxx links to CocHighlightText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 498
CocSymbolDefault xxx links to MoreMsg
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 500
CocPumSearch   xxx links to CocSearch
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 502
CocPumDetail   xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 503
CocPumMenu     xxx links to CocFloating
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 504
CocFloating    xxx ctermbg=253 guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 525
CocPumShortcut xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 505
CocPumDeprecated xxx links to CocStrikeThrough
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 506
CocVirtualText xxx ctermfg=12 guifg=#504945
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 507
CocPumVirtualText xxx links to CocVirtualText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 508
CocInputBoxVirtualText xxx links to CocVirtualText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 509
CocFloatDividingLine xxx links to CocVirtualText
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 510
CocInlineVirtualText xxx ctermfg=244 guifg=#808080
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 512
CocInlineAnnotation xxx links to MoreMsg
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 516
CocListBlackBlack xxx guifg=#282828 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackBlue xxx guifg=#282828 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackGreen xxx guifg=#282828 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackGrey xxx guifg=#282828 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackWhite xxx guifg=#282828 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackCyan xxx guifg=#282828 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackYellow xxx guifg=#282828 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackMagenta xxx guifg=#282828 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlackRed xxx guifg=#282828 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgBlack xxx ctermfg=0 guifg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgBlack xxx ctermbg=0 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListBlueBlack xxx guifg=#458588 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueBlue xxx guifg=#458588 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueGreen xxx guifg=#458588 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueGrey xxx guifg=#458588 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueWhite xxx guifg=#458588 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueCyan xxx guifg=#458588 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueYellow xxx guifg=#458588 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueMagenta xxx guifg=#458588 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListBlueRed xxx guifg=#458588 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgBlue  xxx ctermfg=12 guifg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgBlue  xxx ctermbg=12 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListGreenBlack xxx guifg=#98971a guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenBlue xxx guifg=#98971a guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenGreen xxx guifg=#98971a guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenGrey xxx guifg=#98971a guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenWhite xxx guifg=#98971a guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenCyan xxx guifg=#98971a guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenYellow xxx guifg=#98971a guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenMagenta xxx guifg=#98971a guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreenRed xxx guifg=#98971a guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgGreen xxx ctermfg=10 guifg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgGreen xxx ctermbg=10 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListGreyBlack xxx guifg=#928374 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyBlue xxx guifg=#928374 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyGreen xxx guifg=#928374 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyGrey xxx guifg=#928374 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyWhite xxx guifg=#928374 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyCyan xxx guifg=#928374 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyYellow xxx guifg=#928374 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyMagenta xxx guifg=#928374 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListGreyRed xxx guifg=#928374 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgGrey  xxx ctermfg=248 guifg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgGrey  xxx ctermbg=248 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListWhiteBlack xxx guifg=#a89984 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteBlue xxx guifg=#a89984 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteGreen xxx guifg=#a89984 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteGrey xxx guifg=#a89984 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteWhite xxx guifg=#a89984 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteCyan xxx guifg=#a89984 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteYellow xxx guifg=#a89984 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteMagenta xxx guifg=#a89984 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListWhiteRed xxx guifg=#a89984 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgWhite xxx ctermfg=15 guifg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgWhite xxx ctermbg=15 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListCyanBlack xxx guifg=#689d6a guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanBlue xxx guifg=#689d6a guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanGreen xxx guifg=#689d6a guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanGrey xxx guifg=#689d6a guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanWhite xxx guifg=#689d6a guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanCyan xxx guifg=#689d6a guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanYellow xxx guifg=#689d6a guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanMagenta xxx guifg=#689d6a guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListCyanRed xxx guifg=#689d6a guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgCyan  xxx ctermfg=14 guifg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgCyan  xxx ctermbg=14 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListYellowBlack xxx guifg=#d79921 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowBlue xxx guifg=#d79921 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowGreen xxx guifg=#d79921 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowGrey xxx guifg=#d79921 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowWhite xxx guifg=#d79921 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowCyan xxx guifg=#d79921 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowYellow xxx guifg=#d79921 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowMagenta xxx guifg=#d79921 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListYellowRed xxx guifg=#d79921 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgYellow xxx ctermfg=11 guifg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgYellow xxx ctermbg=11 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListMagentaBlack xxx guifg=#b16286 guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaBlue xxx guifg=#b16286 guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaGreen xxx guifg=#b16286 guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaGrey xxx guifg=#b16286 guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaWhite xxx guifg=#b16286 guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaCyan xxx guifg=#b16286 guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaYellow xxx guifg=#b16286 guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaMagenta xxx guifg=#b16286 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListMagentaRed xxx guifg=#b16286 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgMagenta xxx ctermfg=13 guifg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgMagenta xxx ctermbg=13 guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocListRedBlack xxx guifg=#cc241d guibg=#282828
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedBlue xxx guifg=#cc241d guibg=#458588
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedGreen xxx guifg=#cc241d guibg=#98971a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedGrey xxx guifg=#cc241d guibg=#928374
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedWhite xxx guifg=#cc241d guibg=#a89984
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedCyan xxx guifg=#cc241d guibg=#689d6a
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedYellow xxx guifg=#cc241d guibg=#d79921
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedMagenta xxx guifg=#cc241d guibg=#b16286
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListRedRed  xxx guifg=#cc241d guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 223
CocListFgRed   xxx ctermfg=9 guifg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 225
CocListBgRed   xxx ctermbg=9 guibg=#cc241d
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 226
CocMenuSel     xxx ctermbg=250 guibg=#c8c8c8
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 526
CocListLine    xxx ctermbg=254 guibg=#eaeaea
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 542
CocFloatThumb  xxx ctermbg=249 guibg=#b7b7b7
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 527
CocFloatSbar   xxx links to CocFloating
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 528
CocFloatBorder xxx links to CocFloating
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 539
CocErrorHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocErrorSign   xxx ctermfg=9 guifg=#ff0000
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocErrorVirtualText xxx ctermfg=9 guifg=#ff0000
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocWarningHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocWarningSign xxx ctermfg=130 guifg=#ff922b
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocWarningVirtualText xxx ctermfg=130 guifg=#ff922b
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocInfoHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocInfoSign    xxx ctermfg=11 guifg=#fab005
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocInfoVirtualText xxx ctermfg=11 guifg=#fab005
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocHintHighlight xxx links to CocUnderline
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 562
CocHintSign    xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 567
CocHintVirtualText xxx ctermfg=12 guifg=#15aabf
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocHintFloat   xxx ctermfg=12 ctermbg=253 guifg=#15aabf guibg=#e0e0e0
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocInlayHint   xxx ctermfg=12 ctermbg=248 guifg=#15aabf guibg=Grey
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 236
CocInlayHintParameter xxx links to CocInlayHint
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 583
CocInlayHintType xxx links to CocInlayHint
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 583
CocSemTypeMacro xxx links to Define
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeEnumMember xxx links to Constant
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeComment xxx links to Comment
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeOperator xxx links to Operator
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeProperty xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeClass xxx links to Special
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeStruct xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeRegexp xxx links to String
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeBoolean xxx links to Boolean
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeMethod xxx links to Function
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeNamespace xxx links to Include
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeModifier xxx links to StorageClass
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeNumber xxx links to Number
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeTypeParameter xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeKeyword xxx links to Keyword
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemModDeprecated xxx links to CocDeprecatedHighlight
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeFunction xxx links to Function
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeDecorator xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeEnum xxx links to Type
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeParameter xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeType xxx links to Type
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeString xxx links to String
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeVariable xxx links to Identifier
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeInterface xxx links to Type
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSemTypeEvent xxx links to Keyword
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 617
CocSymbolUnit  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolNumber xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolFunction xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolKey   xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolKeyword xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolReference xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolFolder xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolVariable xxx ctermfg=5 guifg=#6a5acd
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolNull  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolValue xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolConstant xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolText  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolModule xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolPackage xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolClass xxx ctermfg=5 guifg=#6a5acd
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolOperator xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolStruct xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolObject xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolMethod xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolArray xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolEnum  xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolField xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolInterface xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolProperty xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolColor xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolFile  xxx ctermfg=130 guifg=Brown
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolEvent xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolTypeParameter xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolConstructor xxx ctermfg=5 guifg=#6a5acd
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolSnippet xxx ctermfg=2 guifg=SeaGreen
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolBoolean xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolNamespace xxx ctermfg=5 guifg=#6a0dad
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolString xxx ctermfg=1 guifg=Magenta
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSymbolEnumMember xxx ctermfg=6 guifg=DarkCyan
	Last set from ~/.vim/plugged/coc.nvim/plugin/coc.vim line 659
CocSelectedLine xxx cleared
chdir(./RunConfig)
fchdir() to previous dir
chdir(./RunConfig)
fchdir() to previous dir
chdir(/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc)
fchdir() to previous dir
line 10: sourcing "/home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/ui.vim"
finished sourcing /home/upc_gcp/.vim/plugged/coc.nvim/autoload/coc/ui.vim
continuing in coc#api#Call_function
Executing ModeChanged Autocommands for "*"
autocommand call s:Autocmd('ModeChanged', v:event)


chdir(./RunConfig)
fchdir() to previous dir
Executing BufWinLeave Autocommands for "*"
autocommand call s:Autocmd('BufWinLeave', +expand('<abuf>'), bufwinid(+expand('<abuf>')))

chdir(./RunConfig)
fchdir() to previous dir
Executing BufUnload Autocommands for "*"
autocommand call s:Autocmd('BufUnload', +expand('<abuf>'))

Executing VimLeavePre Autocommands for "*"
autocommand call s:VimLeavePre()

Writing viminfo file "/home/upc_gcp/.viminfo"