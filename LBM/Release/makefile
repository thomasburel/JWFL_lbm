################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/User/subdir.mk
-include src/SpecialTreatment/subdir.mk
-include src/Parallelism/subdir.mk
-include src/Mesh/SingleBlock/Patch/subdir.mk
-include src/Mesh/SingleBlock/subdir.mk
-include src/Mesh/MultiBlock/subdir.mk
-include src/IO/subdir.mk
-include src/Core/subdir.mk
-include src/Algorithm/Tools/subdir.mk
-include src/Algorithm/LowOrder/Boundaries/subdir.mk
-include src/Algorithm/LowOrder/subdir.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: LBM

# Tool invocations
LBM: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	mpic++ -o "LBM" $(OBJS) $(USER_OBJS) $(LIBS) -lcgns -lhdf5 -lboost_serialization
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) LBM
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
