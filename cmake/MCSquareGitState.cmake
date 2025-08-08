execute_process(COMMAND git log --pretty=format:'%h' -n 1
				WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE MCsquare_GIT_REV)
if (MCsquare_GIT_REV)
    execute_process(
        COMMAND git describe --exact-match --tags
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE MCsquare_GIT_TAG ERROR_QUIET)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE MCsquare_GIT_BRANCH)

    string(STRIP "${MCsquare_GIT_REV}" MCsquare_GIT_REV)
    string(SUBSTRING "${MCsquare_GIT_REV}" 1 7 MCsquare_GIT_REV)
    string(STRIP "${MCsquare_GIT_TAG}" MCsquare_GIT_TAG)
    string(STRIP "${MCsquare_GIT_BRANCH}" MCsquare_GIT_BRANCH)
endif()

