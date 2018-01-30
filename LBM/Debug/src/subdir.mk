################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Doxygen.cpp \
../src/LBM.cpp 

OBJS += \
./src/Doxygen.o \
./src/LBM.o 

CPP_DEPS += \
./src/Doxygen.d \
./src/LBM.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -o "$@" "$<" -O0 -g3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)"
	@echo 'Finished building: $<'
	@echo ' '


