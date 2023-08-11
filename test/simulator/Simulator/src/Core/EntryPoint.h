#pragma once

#ifdef PLATFORM_WINDOWS
extern temp::Application* temp::Creat_application();
int main(int argc, char** argv) {
	auto app = temp::Creat_application();
	app->Run();
	delete app;
}

#endif