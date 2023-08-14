#pragma once

#ifdef PLATFORM_WINDOWS
	#ifdef BUILD_DLL
		#define C_API __declspec(dllexport)
	#else
		#define C_API __declspec(dllimport)
	#endif // BUILE_DLL
#else
	#error not work on this!
#endif // PLATFORM_WINDOWS
