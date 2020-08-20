print("Loading libraries")
import asyncio
import os
import sys
import time

import aiohttp
import jwt
from gidgethub.aiohttp import GitHubAPI
from github3 import login


def get_jwt(app_id):

    pem_file = os.getenv("PEM_FILE")
    if pem_file is None:
        raise Exception("PEM_FILE not specified")

    payload = {
        "iat": int(time.time()),
        "exp": int(time.time()) + (10 * 60),
        "iss": app_id,
    }
    encoded = jwt.encode(payload, pem_file, algorithm="RS256")
    bearer_token = encoded.decode("utf-8")

    return bearer_token


async def get_installation(gh, jwt, username):
    async for installation in gh.getiter(
        "/app/installations",
        jwt=jwt,
        accept="application/vnd.github.machine-man-preview+json",
    ):
        if installation["account"]["login"] == username:
            return installation

    raise ValueError(f"Can't find installation by that user: {username}")


async def get_installation_access_token(gh, jwt, installation_id):
    # doc: https: // developer.github.com/v3/apps/#create-a-new-installation-token

    access_token_url = (
        f"https://api.github.com/app/installations/{installation_id}/access_tokens"
    )
    response = await gh.post(
        access_token_url,
        data=b"",
        jwt=jwt,
        accept="application/vnd.github.machine-man-preview+json",
    )
    # example response
    # {
    #   "token": "v1.1f699f1069f60xxx",
    #   "expires_at": "2016-07-11T22:14:10Z"
    # }

    return response


async def main(release):
    print("Determining version")
    version = open("dist/version").read().strip()
    print("Version:", version)
    tarball="truchas-%s-Linux.tar.bz2" % version
    print("Tarball:", tarball)
    print("Authenticating")
    async with aiohttp.ClientSession() as session:
        app_id = os.getenv("GH_APP_ID")
        if app_id is None:
            print("GH_APP_ID not specified, skipping")
            return

        jwt = get_jwt(app_id)
        gh = GitHubAPI(session, "truchas")

        print("Obtaining installation")
        try:
            installation = await get_installation(gh, jwt, "truchas")

        except ValueError as ve:
            # Raised if the user did not install the GitHub App
            print(ve)
        else:
            print("Obtaining installation token")
            access_token = await get_installation_access_token(
                gh, jwt=jwt, installation_id=installation["id"]
            )

            print("Token obtained")
            gh = login("TruchasUploader", access_token["token"])
            repo = gh.repository("truchas", "truchas_releases")


            print("Release information:")

            url_base = "https://github.com/truchas/truchas_releases/releases/download"
            tarball_url = "{url_base}/{version}/{tarball}".format(
                    url_base=url_base,
                    version=version,
                    tarball=tarball,
                    )
            release_url="https://github.com/truchas/truchas_releases/releases/tag/{version}".format(version=version)
            tarball_path = "dist/" + tarball
            tarball_size = len(open(tarball_path, "rb").read())
            release_name = "Release version %s" % version
            release_body="""\
Truchas binary release. Version: {version}

| Filename | Size |
| -------- | ---- |
| [{tarball}]({tarball_url}) | {t_size_MB} MB |

""".format(
        tarball=tarball,
        tarball_url=tarball_url,
        version=version,
        t_size_MB="%.1f" % (tarball_size / 1024.**2),
    )

            print("Release name:", release_name)
            print("Release body:")
            print(release_body)
            print("Release url:")
            print(release_url)
            print("Truchas tarball url:")
            print(tarball_url)

            if not release:
                print("Not uploading release.")
                return

            print("Creating a release")
            r = repo.create_release(version,
                    name=release_name,
                    body=release_body,
                    draft=False)
            print("Release URL:")
            print(r.html_url)
            print("Uploading binay tarball")
            f = open(tarball_path, "rb")
            s = r.upload_asset("application/x-bzip2", tarball, f)
            print("Uploaded:")
            print(s.browser_download_url)



release = False
if (len(sys.argv) == 2 and sys.argv[1] == "master"):
    # If the first argument is "master", it is a development version, upload it
    release = True
if (len(sys.argv) == 3 and sys.argv[2] != ""):
    # If the second argument is a tag, it is a release version, upload it
    release = True
print("Upload binary:", release)

asyncio.run(main(release))
