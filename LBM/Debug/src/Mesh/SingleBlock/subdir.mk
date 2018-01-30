################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Mesh/SingleBlock/Block.cpp \
../src/Mesh/SingleBlock/Block2D.cpp \
../src/Mesh/SingleBlock/Cell2D.cpp \
../src/Mesh/SingleBlock/Node2D.cpp \
../src/Mesh/SingleBlock/NodeArrays.cpp 

OBJS += \
./src/Mesh/SingleBlock/Block.o \
./src/Mesh/SingleBlock/Block2D.o \
./src/Mesh/SingleBlock/Cell2D.o \
./src/Mesh/SingleBlock/Node2D.o \
./src/Mesh/SingleBlock/NodeArrays.o 

CPP_DEPS += \
./src/Mesh/SingleBlock/Block.d \
./src/Mesh/SingleBlock/Block2D.d \
./src/Mesh/SingleBlock/Cell2D.d \
./src/Mesh/SingleBlock/Node2D.d \
./src/Mesh/SingleBlock/NodeArrays.d 


# Each subdirectory must supply rules for building sources it contributes
src/Mesh/SingleBlock/%.o: ../src/Mesh/SingleBlock/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -o "$@" "$<" -O0 -g3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)"
	@echo 'Finished building: $<'
	@echo ' '


