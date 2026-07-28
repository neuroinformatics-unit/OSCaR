# PyRAT

[PyRAT](https://www.scionics.com/pyrat.html) is a popular animal facility management software.

To access data in a PyRAT instance via OSCaR, there are some initial setup steps:

## Client tokens

This step should be done by an administrator of your institution's PyRAT instance, and only needs to be done once.

You will need to contact Scionics (the creators of PyRAT) to approve access for OSCaR:

- Ask them to setup a new API client called 'OSCaR' on your pyRAT server; they should send you a 'client token'.
- The new API client will appear under `ADMINISTRATION > API` in pyRAT, and an administrator will need to enable it from there.

## User tokens

Once the client is enabled, any user can request a user token by logging into PyRAT, going to `ADMINISTRATION > API` and clicking the 'Request access' button on the `OSCaR` row.

## Providing tokens to OSCaR

Once all tokens have been created (as above), they can be provided to OSCaR by setting the following environment variables:

- `PYRAT_URL`: the url of the pyRAT instance you want to access
- `PYRAT_CLIENT_TOKEN`: the client token sent to you by Scionics
- `PYRAT_USER_TOKEN`: the user token created via 'Request access' in PyRAT

Set these variables however you prefer - just make sure they are _never_ checked in to version control (tokens should be kept secret). One convenient method is to create a `.env` file and use [`python-dotenv`](https://github.com/theskumar/python-dotenv) to read them. We'll show an example of that below.

## Pulling data from the PyRAT api

## Standardising data from the PyRAT api
