ORDER 	?= 4095

SRCS    := main.c gbst.c gbstlock.c
TARGET  := GreenBST
TREE	:= -DGBST -D__PREALLOCGNODES=${ORDER}


include  ../common/common.mk

withorder: clean ${TARGET}.energy
	@cp ${TARGET}.energy ${TARGET}.energy.${ORDER}
