#!/bin/sh

# set install path
HOMEDIR=~
PREFIX="$HOMEDIR/bin/"
VIMDIR="$HOMEDIR/.vim/"

# get relevant scripts
echo "Getting macvim-load-line from googlecode repository"
curl -s http://agpy.googlecode.com/svn/trunk/macvim-skim/macvim-load-line > ${PREFIX}/macvim-load-line
chmod +x ${PREFIX}/macvim-load-line
echo "Getting WhichTab.vim from googlecode repository"
curl -s http://agpy.googlecode.com/svn/trunk/macvim-skim/WhichTab.vim > ${VIMDIR}/plugin/WhichTab.vim

if [ ! -d /Applications/Skim.app ]
then
    # get Skim.app 
    # should work, but doesn't curl -L http://sourceforge.net/projects/skim-app/files/latest/download?source=files -o Skim.dmg
    echo "Did not find /Applications/Skim.app.  Downloading and installing 1.3.19."
    curl -L "http://downloads.sourceforge.net/project/skim-app/Skim/Skim-1.3.19/Skim-1.3.19.dmg?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fskim-app%2F&ts=1326916683&use_mirror=superb-dca2" -o Skim.dmg
    hdid Skim.dmg
    cp -r /Volumes/Skim/Skim.app /Applications
    hdiutil eject /Volumes/Skim
fi

# Add lines to .vimrc 
if [[ `grep 'map ,[rvpmt]' ${HOMEDIR}/.vimrc` == "" ]]; 
then
    echo "Did not find mappings in ${HOMEDIR}/.vimrc.  Appending them."
    echo \" Activate skim >> ${HOMEDIR}/.vimrc
    echo 'map ,v :w<CR>:silent !/Applications/Skim.app/Contents/SharedSupport/displayline -r <C-r>=line(".")<CR> %<.pdf %<CR><CR>' >> ${HOMEDIR}/.vimrc
    echo 'map ,p :w<CR>:silent !pdflatex -synctex=1 --interaction=nonstopmode %:p <CR>:silent !/Applications/Skim.app/Contents/SharedSupport/displayline -r <C-r>=line(".")<CR> %<.pdf %<CR><CR>' >> ${HOMEDIR}/.vimrc
    echo 'map ,m :w<CR>:silent !make <CR>:silent !/Applications/Skim.app/Contents/SharedSupport/displayline -r <C-r>=line(".")<CR> %<.pdf %<CR><CR>' >> ${HOMEDIR}/.vimrc
    echo \" Reactivate VIM >> ${HOMEDIR}/.vimrc
    echo 'map ,r :w<CR>:silent !/Applications/Skim.app/Contents/SharedSupport/displayline -r <C-r>=line(".")<CR> %<.pdf %<CR>:silent !osascript -e "tell application \"MacVim\" to activate" <CR><CR>' >> ${HOMEDIR}/.vimrc
    echo 'map ,t :w<CR>:silent !pdflatex -synctex=1 --interaction=nonstopmode %:p <CR>:silent !/Applications/Skim.app/Contents/SharedSupport/displayline -r <C-r>=line(".")<CR> %<.pdf %<CR>:silent !osascript -e "tell application \"MacVim\" to activate" <CR><CR>' >> ${HOMEDIR}/.vimrc
fi

# set Skim settings
defaults write net.sourceforge.skim-app.skim SKTeXEditorCommand "macvim-load-line"
defaults write net.sourceforge.skim-app.skim SKTeXEditorArguments "'\"%file\" %line'"
