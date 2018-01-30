################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/SpecialTreatment/PorousMedia.cpp 

OBJS += \
./src/SpecialTreatment/PorousMedia.o 

CPP_DEPS += \
./src/SpecialTreatment/PorousMedia.d 


# Each subdirectory must supply rules for building sources it contributes
src/SpecialTreatment/%.o: ../src/SpecialTreatment/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -o "$@" "$<" -O3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)"
	@echo 'Finished building: $<'
	@echo ' '


