################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/User/UserContactAngle.cpp \
../src/User/UserForce.cpp \
../src/User/UserInit.cpp \
../src/User/UserMesh.cpp \
../src/User/UserParameters.cpp \
../src/User/UserPatchBc.cpp 

OBJS += \
./src/User/UserContactAngle.o \
./src/User/UserForce.o \
./src/User/UserInit.o \
./src/User/UserMesh.o \
./src/User/UserParameters.o \
./src/User/UserPatchBc.o 

CPP_DEPS += \
./src/User/UserContactAngle.d \
./src/User/UserForce.d \
./src/User/UserInit.d \
./src/User/UserMesh.d \
./src/User/UserParameters.d \
./src/User/UserPatchBc.d 


# Each subdirectory must supply rules for building sources it contributes
src/User/%.o: ../src/User/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -o "$@" "$<" -O0 -g3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)"
	@echo 'Finished building: $<'
	@echo ' '


