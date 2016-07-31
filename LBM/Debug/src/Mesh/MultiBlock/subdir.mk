################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Mesh/MultiBlock/MultiBlock.cpp \
../src/Mesh/MultiBlock/MultiBlock2D.cpp 

OBJS += \
./src/Mesh/MultiBlock/MultiBlock.o \
./src/Mesh/MultiBlock/MultiBlock2D.o 

CPP_DEPS += \
./src/Mesh/MultiBlock/MultiBlock.d \
./src/Mesh/MultiBlock/MultiBlock2D.d 


# Each subdirectory must supply rules for building sources it contributes
src/Mesh/MultiBlock/%.o: ../src/Mesh/MultiBlock/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -O0 -g3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


