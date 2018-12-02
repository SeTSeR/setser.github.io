---
layout: post
author: Sergey
title: "Get track data from D-Bus"
comments: true
---

An example of getting info about currently playing track from D-Bus using MPRIS. This code was used in [AnnicomScrobbler](https://github.com/SeTSeR/AnnicomScrobbler.git).

```c
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>
#include <glib/gprintf.h>
#include <gio/gio.h>

#include <curl/curl.h>

struct song {
  gchar *title;
  gchar *artist;
  gchar *genre;
};

void get_song(char* address, struct song *currsong) {
  GDBusProxy *proxy;
  GDBusConnection *conn;
  GError *error = NULL;
  GVariant *answer, *value;
  GVariantIter *iter;
  gchar *key;

  conn = g_bus_get_sync(G_BUS_TYPE_SESSION, NULL, &error);
  g_assert_no_error(error);

  proxy = g_dbus_proxy_new_sync(conn,
                                G_DBUS_PROXY_FLAGS_NONE,
                                NULL,
                                address,
                                "/org/mpris/MediaPlayer2",
                                "org.freedesktop.DBus.Properties",
                                NULL,
                                &error);
  if(error != NULL) {
    g_object_unref(proxy);
    g_object_unref(conn);
    return;
  }

  answer = g_dbus_proxy_call_sync(proxy,
                                  "Get",
                                  g_variant_new("(ss)",
                                                "org.mpris.MediaPlayer2.Player",
                                                "Metadata"),
                                  G_DBUS_CALL_FLAGS_NONE,
                                  -1,
                                  NULL,
                                  &error);
  if((error != NULL) || (answer == NULL)) {
    if(answer != NULL) {
      g_variant_unref(answer);
    }
    g_object_unref(proxy);
    g_object_unref(conn);
    return;
  }

  GVariant *tmp;
  g_variant_get(answer,"(v)", &tmp);
  g_variant_get(tmp, "a{sv}", &iter);
  while(g_variant_iter_loop(iter, "{sv}", &key, &value)) {
    if(strncmp(key, "xesam:title", 11) == 0) {
      g_variant_get(value, "s", &(currsong->title));
    }
    else if(strncmp(key, "xesam:genre", 11) == 0) {
      if(g_variant_is_container(value)) {
        g_variant_get_child(value, 0, "s", &(currsong->genre));
      }
      else {
        g_variant_get(value, "s", &(currsong->genre));
      }
    }
    else if(strncmp(key, "xesam:artist", 12) == 0) {
      if(g_variant_is_container(value)) {
        g_variant_get_child(value, 0, "s", &(currsong->artist));
      }
      else {
        g_variant_get(value, "s", &(currsong->artist));
      }
    }
  }

  g_variant_unref(tmp);
  g_variant_unref(answer);

  g_object_unref(proxy);
  g_object_unref(conn);
}
```

Usage:

```c
get_song(currsong);
printf("song title: %s\ngenre: %s, artist: %s", currsong->title, currsong->genre, currsong->artist);
```
