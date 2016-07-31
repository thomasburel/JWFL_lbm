################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/IO/CGNS.cpp \
../src/IO/TecplotIO.cpp \
../src/IO/WriterManager.cpp 

OBJS += \
./src/IO/CGNS.o \
./src/IO/TecplotIO.o \
./src/IO/WriterManager.o 

CPP_DEPS += \
./src/IO/CGNS.d \
./src/IO/TecplotIO.d \
./src/IO/WriterManager.d 


# Each subdirectory must supply rules for building sources it contributes
src/IO/%.o: ../src/IO/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -O0 -g3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


