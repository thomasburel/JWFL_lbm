################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Parallelism/MpiManager.cpp \
../src/Parallelism/ParallelManager.cpp 

OBJS += \
./src/Parallelism/MpiManager.o \
./src/Parallelism/ParallelManager.o 

CPP_DEPS += \
./src/Parallelism/MpiManager.d \
./src/Parallelism/ParallelManager.d 


# Each subdirectory must supply rules for building sources it contributes
src/Parallelism/%.o: ../src/Parallelism/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -O3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


