#pragma once

#include "core.h"
namespace temp {

	class C_API Application
	{
	public:
		Application();
		virtual ~Application();	
		void Run();
	};
	Application* Creat_application();

}