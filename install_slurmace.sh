#!/bin/sh
curl -O https://raw.githubusercontent.com/asqwerty666/acenip/main/lib/SLURMACE.pm
mkdir -p $HOME/.local/lib/perl5
mv SLURMACE.pm $HOME/.local/lib/perl5/
export PERL5LIB=$PERL5LIB:$HOME/.local/lib/perl5/
echo 'export PERL5LIB=$PERL5LIB:$HOME/.local/lib/perl5/' >> $HOME/.bash_profile
