#!/bin/bash

. ./edafunctions.sh --source-only

export AXIS="z"

export SEEDVAL=1.3
export SYSYEMNAME="AuHg+"
posteda

export SEEDVAL=1.3
export SYSYEMNAME="AuCn+"
posteda

export SEEDVAL=1.3
export SYSYEMNAME="AuPb+"
posteda

export SEEDVAL=1.3
export SYSYEMNAME="AuFl+"
posteda

export SEEDVAL=1.4
export SYSYEMNAME="AuRn+"
posteda

export SEEDVAL=1.4
export SYSYEMNAME="AuOg+"
posteda
