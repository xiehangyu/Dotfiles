# Enable Powerlevel10k instant prompt. Should stay close to the top of ~/.zshrc.
# Initialization code that may require console input (password prompts, [y/n]
# confirmations, etc.) must go above this block, everything else may go below.
if [[ -r "${XDG_CACHE_HOME:-$HOME/.cache}/p10k-instant-prompt-${(%):-%n}.zsh" ]]; then
  source "${XDG_CACHE_HOME:-$HOME/.cache}/p10k-instant-prompt-${(%):-%n}.zsh"
fi

#!/bin/zsh
# If you come from bash you might have to change your $PATH.
 export PATH=$HOME/bin:/usr/local/bin:$PATH

# Path to your oh-my-zsh installation.
 export ZSH=/home/petergu/.oh-my-zsh

# Set name of the theme to load. Optionally, if you set this to "random"
# it'll load a random theme each time that oh-my-zsh is loaded.
# See https://github.com/robbyrussell/oh-my-zsh/wiki/Themes
if [[ $TERM == *256color ]];then
	#ZSH_THEME="sobole"
	#murilasso
	#simonoff
	#ZSH_THEME="robbyrussell"
	ZSH_THEME="powerlevel10k/powerlevel10k"
	#ZSH_THEME="agnoster"
	#ZSH_THEME="random"
else
	ZSH_THEME="robbyrussell"
	#ZSH_THEME="sobole"
fi
# Set list of themes to load
# Setting this variable when ZSH_THEME=random
# cause zsh load theme from this variable instead of
# looking in ~/.oh-my-zsh/themes/
# An empty array have no effect
# ZSH_THEME_RANDOM_CANDIDATES=( "robbyrussell" "agnoster" )

# Uncomment the following line to use case-sensitive completion.
# CASE_SENSITIVE="true"

# Uncomment the following line to use hyphen-insensitive completion. Case
# sensitive completion must be off. _ and - will be interchangeable.
# HYPHEN_INSENSITIVE="true"

# Uncomment the following line to disable bi-weekly auto-update checks.
 DISABLE_AUTO_UPDATE="true"

# Uncomment the following line to change how often to auto-update (in days).
# export UPDATE_ZSH_DAYS=13

# Uncomment the following line to disable colors in ls.
# DISABLE_LS_COLORS="true"

# Uncomment the following line to disable auto-setting terminal title.
# DISABLE_AUTO_TITLE="true"

# Uncomment the following line to enable command auto-correction.
 ENABLE_CORRECTION="true"

# Uncomment the following line to display red dots whilst waiting for completion.
 COMPLETION_WAITING_DOTS="true"

# Uncomment the following line if you want to disable marking untracked files
# under VCS as dirty. This makes repository status check for large repositories
# much, much faster.
# DISABLE_UNTRACKED_FILES_DIRTY="true"

# Uncomment the following line if you want to change the command execution time
# stamp shown in the history command output.
# The optional three formats: "mm/dd/yyyy"|"dd.mm.yyyy"|"yyyy-mm-dd"
 HIST_STAMPS="mm/dd/yyyy"

# Would you like to use another custom folder than $ZSH/custom?
# ZSH_CUSTOM=/path/to/new-custom-folder

# Which plugins would you like to load? (plugins can be found in ~/.oh-my-zsh/plugins/*)
# Custom plugins may be added to ~/.oh-my-zsh/custom/plugins/
# Example format: plugins=(rails git textmate ruby lighthouse)
# Add wisely, as too many plugins slow down shell startup.
plugins=(
git
#command-not-found
#zsh-autosuggestions
zsh-syntax-highlighting
#lighthouse
)


#[[ $TERM == 'xterm-256color' ]] && source $ZSH/oh-my-zsh.sh
source $ZSH/oh-my-zsh.sh

# User configuration

# export MANPATH="/usr/local/man:$MANPATH"

# You may need to manually set your language environment
# export LANG=en_US.UTF-8

# Preferred editor for local and remote sessions
# if [[ -n $SSH_CONNECTION ]]; then
#   export EDITOR='vim'
# else
#   export EDITOR='mvim'
# fi

# Compilation flags
# export ARCHFLAGS="-arch x86_64"

export DEFAULT_USER="petergu"

# ssh
# export SSH_KEY_PATH="~/.ssh/rsa_id"

# Set personal aliases, overriding those provided by oh-my-zsh libs,
# plugins, and themes. Aliases can be placed here, though oh-my-zsh
# users are encouraged to define aliases within the ZSH_CUSTOM folder.
# For a full list of active aliases, run `alias`.
#
# Example aliases
# alias zshconfig="mate ~/.zshrc"
# alias ohmyzsh="mate ~/.oh-my-zsh"
alias vi="vim --noplugin"
#alias nvi="nvim --noplugin"
alias ll="ls -alh"
#alias lll="exa -abghHliS"

#give up save/load and enjoy real game.
#alias nha="~/Python/nha/main.py"
#source /home/petergu/Downloads/zsh-syntax-highlighting/zsh-syntax-highlighting.zsh

#function s(){
    #local spath="/usr/local/bin/s"
    #if [[ ${1:0:1} != '-' ]]
    #then
        #$spath $* | less 
    #else
        #$spath $*
    #fi
#}
#alias s='echo -ne '\n'|s'

#alias matlab='/usr/local/MATLAB/R2015b/bin/matlab'
#alias matcmd='matlab -nodisplay -nosplash -nojvm'
#alias zh='source ~petergu/Widgets/zh.sh'
alias pg='watch -n 1 progress'
#some legacy from my old good slow hdd mbp8,1
#alias note='nvim ~petergu/Widgets/news.txt'
#function note(){
	#local notepath="/home/petergu/Widgets"
	#local notefile="news.txt"
	#local edit="nvim"
	#mv -f "$notepath/$notefile" "$notepath/.$notefile.tmp"
	#if [[ "$1" == "d" ]]
	#then
		#date > "$notepath/$notefile"
	#fi
	#cat "$notepath/.$notefile.tmp" >> "$notepath/$notefile"
	##rm "$notepath/.$notefile.tmp"
	#$edit "$notepath/$notefile"
	
#}
#alias dnote='note d'
#alias origin='wine "C:\Program Files (x86)\OriginLab\Origin2017\Origin94.exe"'
#export WINEDEBUG=-all
export LC_TIME="en_US.UTF-8"
#export XIM="fcitx"
#export XIM_PROGRAM="fcitx"
#export XMODIFIERS="@im=fcitx"
#export GTK_IM_MODULE="fcitx"
#export QT_IM_MODULE="fcitx"
#export http_proxy=http://xxx:xxx@192.168.2.49:8080
alias vrc="nvim ~/.vimrc"
alias zrc="nvim ~/.zshrc"
#alias ai="sudo apt install"
#alias arm="sudo apt remove"
#alias apg="sudo apt purge"
#alias di="sudo dpkg -i"
#alias drm="sudo dpkg -r"
#alias dpg="sudo dpkg -P"
#alias af="sudo apt -f install"
#alias afi="sudo apt-fast install"
#alias au="sudo apt update"
alias pi="sudo pacman -S"
alias prm="sudo pacman -R"
alias ppi="yay -S"
alias py3="python3"
alias py2="python2"
alias adown="axel -n 10"
alias batt="sudo tlp batt"

alias autophy="/home/petergu/PhysicsExp/Core/autophy"

alias news="vim /home/petergu/Widgets/news.txt"

#a better script: ~/Widget/deepin.sh
#alias dwps="sudo xchroot /mnt 'su petergu -c wps'"
#alias dwpp="sudo xchroot /mnt 'su petergu -c wpp'"
#alias det="sudo xchroot /mnt 'su petergu -c et'"
#alias dtim="sudo xchroot /mnt 'su petergu -c \"env WINEPREFIX=/home/petergu/.deepinwine/Deepin-TIM/ deepin-wine wineboot\"'"
#alias dwechat="sudo xchroot /mnt 'su petergu -c \"env WINEPREFIX=/home/petergu/.deepinwine/Deepin-WeChat/ deepin-wine wineboot\"'"

#export TEXMACS_PATH='/home/petergu/src/TeXmacs-1.99.7-11201-x86_64-pc-linux-gnu/TeXmacs'
#export PATH=$TEXMACS_PATH/bin:$PATH

#export GOROOT=/usr/local/go
#export PATH=$PATH:$GOROOT/bin

#export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH

export PATH=$PATH:~/.gem/ruby/2.6.0/bin

function light()
{
	export http_proxy="http://user:pass@example.com:port/"; export https_proxy=$http_proxy;
	exec zsh
}
function ssr()
{
	export http_proxy="http://127.0.0.1:12333"
	export https_proxy="http://127.0.0.1:12333"
	exec zsh
}
function v2()
{
	export http_proxy="http://127.0.0.1:7890"
	export https_proxy="http://127.0.0.1:7890"
	exec zsh
}

[ -f ~/.fzf.zsh ] && source ~/.fzf.zsh
eval `dircolors ~/.dircolors`
# to avoid confusion between wsl and VM
# export PS1='${ret_status} SSH%{$fg[cyan]%}:%c%{$reset_color%} $(git_prompt_info)'

# LLVM
#export PATH=/run/media/petergu/ssid/petergu/llvm-install/bin:$PATH
export PATH=/home/petergu/MyHome/src/llvm-install/bin:$PATH
export PATH=/opt/android-sdk/platform-tools:$PATH
#export PATH=/run/media/petergu/ssid/petergu/riscv/build/bin:$PATH
#export PATH=/opt/riscv/bin:$PATH
#export PATH=/opt/riscv/riscv64-unknown-linux-gnu/bin/:$PATH


eval $(thefuck --alias)

# To customize prompt, run `p10k configure` or edit ~/.p10k.zsh.
[[ ! -f ~/.p10k.zsh ]] || source ~/.p10k.zsh
