"""SimForge ↔ VMD bridge protocol, version ``simforge-vmd/1``.

Line-oriented ASCII over loopback TCP; one message per ``\\n``-terminated line;
at most :data:`MAX_LINE` bytes; fields separated by single spaces.  Parsing is
allow-listed and deterministic: nothing received is ever executed.

VMD bridge → SimForge::

    HELLO simforge-vmd/1 <token>
    LOADING
    READY <molid> <numatoms> <numframes> <vmd_version>
    FRAME <frame> <ack_sequence|NONE>
    ERROR <code> [message…]
    BYE

SimForge → VMD bridge::

    WELCOME simforge-vmd/1
    GOTO <sequence> <frame>          # sequence 0 = adapter-internal initial sync
    BYE

Frames are VMD frame indices after the structure file's own coordinate frame
was removed, i.e. exactly the zero-based display frames of the ReviewDataset.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

PROTOCOL = "simforge-vmd/1"
MAX_LINE = 512
INITIAL_SYNC_SEQUENCE = 0

_INT = re.compile(r"^(0|[1-9][0-9]{0,9})$")
_TOKEN = re.compile(r"^[0-9a-f]{32}$")
_CODE = re.compile(r"^[a-z_]{1,40}$")
_VERSION = re.compile(r"^[A-Za-z0-9._+-]{1,40}$")
_PRINTABLE = re.compile(r"^[\x20-\x7e]*$")


class ProtocolError(ValueError):
    """A line that is not a valid message of this protocol version."""


@dataclass(frozen=True)
class Message:
    kind: str
    fields: tuple = ()


def _int(text: str, what: str) -> int:
    if not _INT.match(text):
        raise ProtocolError(f"{what} must be a non-negative integer, got {text[:20]!r}")
    return int(text)


def parse_bridge_line(raw: bytes | str) -> Message:
    """Parse one line received *from* the VMD bridge."""
    if isinstance(raw, bytes):
        if len(raw) > MAX_LINE:
            raise ProtocolError("line too long")
        try:
            raw = raw.decode("ascii")
        except UnicodeDecodeError as exc:
            raise ProtocolError("non-ASCII line") from exc
    line = raw.rstrip("\n").rstrip("\r")
    if len(line) > MAX_LINE:
        raise ProtocolError("line too long")
    if not line or not _PRINTABLE.match(line):
        raise ProtocolError("empty or non-printable line")
    parts = line.split(" ")
    cmd, args = parts[0], parts[1:]
    if cmd == "HELLO":
        if len(args) != 2:
            raise ProtocolError("HELLO needs <protocol> <token>")
        if args[0] != PROTOCOL:
            raise ProtocolError(f"unsupported protocol {args[0][:40]!r} (expected {PROTOCOL})")
        if not _TOKEN.match(args[1]):
            raise ProtocolError("malformed token")
        return Message("HELLO", (args[0], args[1]))
    if cmd in ("LOADING", "BYE"):
        if args:
            raise ProtocolError(f"{cmd} takes no arguments")
        return Message(cmd)
    if cmd == "READY":
        if len(args) != 4:
            raise ProtocolError("READY needs <molid> <numatoms> <numframes> <vmd_version>")
        if not _VERSION.match(args[3]):
            raise ProtocolError("malformed VMD version")
        return Message("READY", (_int(args[0], "molid"), _int(args[1], "numatoms"),
                                 _int(args[2], "numframes"), args[3]))
    if cmd == "FRAME":
        if len(args) != 2:
            raise ProtocolError("FRAME needs <frame> <ack_sequence|NONE>")
        ack = None if args[1] == "NONE" else _int(args[1], "ack_sequence")
        return Message("FRAME", (_int(args[0], "frame"), ack))
    if cmd == "ERROR":
        if not args or not _CODE.match(args[0]):
            raise ProtocolError("ERROR needs a lower-case code")
        return Message("ERROR", (args[0], " ".join(args[1:])[:300]))
    raise ProtocolError(f"unknown message {cmd[:20]!r}")


def welcome() -> bytes:
    return f"WELCOME {PROTOCOL}\n".encode()


def goto(sequence: int, frame: int) -> bytes:
    for v, what in ((sequence, "sequence"), (frame, "frame")):
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            raise ProtocolError(f"{what} must be a non-negative integer")
    return f"GOTO {sequence} {frame}\n".encode()


def bye() -> bytes:
    return b"BYE\n"


def ack_of(msg: Message) -> Optional[int]:
    return msg.fields[1] if msg.kind == "FRAME" else None
