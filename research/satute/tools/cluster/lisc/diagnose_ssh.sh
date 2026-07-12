#!/usr/bin/env bash
set -u

HOSTS=(login01.lisc.univie.ac.at login02.lisc.univie.ac.at lisc.univie.ac.at)

echo "SSH aliases"
ssh -G lisc 2>/dev/null | awk '/^(host|hostname|user|identityfile|identitiesonly|proxycommand|proxyjump) / {print}'
echo
ssh -G lisc02 2>/dev/null | awk '/^(host|hostname|user|identityfile|identitiesonly|proxycommand|proxyjump) / {print}'

echo
echo "Key file"
ls -l ~/.ssh/lisc_sakalli_ed25519 ~/.ssh/lisc_sakalli_ed25519.pub 2>&1
ssh-keygen -lf ~/.ssh/lisc_sakalli_ed25519.pub 2>&1

echo
echo "DNS"
for host in "${HOSTS[@]}" google.com; do
  echo "== ${host}"
  dscacheutil -q host -a name "${host}" 2>&1 || true
done

echo
echo "Port 22"
for host in "${HOSTS[@]}"; do
  echo "== ${host}"
  nc -vz -G 5 "${host}" 22 2>&1 || true
done

echo
echo "Batch SSH"
ssh -o BatchMode=yes -o ConnectTimeout=8 lisc 'hostname; whoami; pwd; id' 2>&1 || true
ssh -o BatchMode=yes -o ConnectTimeout=8 lisc02 'hostname; whoami; pwd; id' 2>&1 || true
