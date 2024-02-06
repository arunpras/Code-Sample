#!/usr/bin/env bash
set -eo pipefail


function notify_app() {
    local MSG="$1"
    local SEVERITY=${2:-MESSAGE}
    if [ -n "$MSG" ]; then
        log_message "Sending a message to the app ($SEVERITY): $MSG"
        if hash send_message 2>/dev/null; then
            send_message --message "$MSG" --severity "$SEVERITY"
        else
            echo "Did not send a message to the app, check if send_message is in PATH env"
        fi
    fi
}

function notify_value_error_from() {
    FILE="${1}"
    MSG=$(grep -e "^ValueError" "$FILE" | sed -e "s/^ValueError: //")
    notify_app "$MSG" "HARD_ERROR"
}
