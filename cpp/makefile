CURRENT_DIR := ${CURDIR}
LOCAL_EIGEN_PATH := /eigen3
EIGEN_PATH :=$(addprefix $(CURDIR),$(LOCAL_EIGEN_PATH))

all: 
	@echo --- Building app ${EIGEN_PATH}
	g++ app.cpp -I${EIGEN_PATH} -o app
