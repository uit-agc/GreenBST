ORDER ?= 64

SRCS    := main.c
TARGET  := CBTree
TREE	:= -DDEFAULT_ORDER=${ORDER} -DCBTREE

include  ../common/common.mk

withorder: clean ${TARGET}.energy
	@cp ${TARGET}.energy ${TARGET}.energy.${ORDER}
