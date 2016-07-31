################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/User/UserForce.cpp \
../src/User/UserInit.cpp \
../src/User/UserMesh.cpp \
../src/User/UserParameters.cpp 

OBJS += \
./src/User/UserForce.o \
./src/User/UserInit.o \
./src/User/UserMesh.o \
./src/User/UserParameters.o 

CPP_DEPS += \
./src/User/UserForce.d \
./src/User/UserInit.d \
./src/User/UserMesh.d \
./src/User/UserParameters.d 


# Each subdirectory must supply rules for building sources it contributes
src/User/%.o: ../src/User/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -O3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


